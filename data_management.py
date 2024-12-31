import pandas as pd
import os
import re

import xyz2tab as x2t
import warnings

# Suppress warnings from a specific file
warnings.filterwarnings("ignore", module="xyz2tab")


# --- This method gets the lines of a file after a specific marker --- #
def read_after_last_instance(filepath, marker):
    with open(filepath, "r") as file:
        lines = file.readlines()

    # Find the last occurrence of the marker
    last_occurrence = None
    for i, line in enumerate(lines):
        if marker in line:
            last_occurrence = i

    # Read lines after the last occurrence
    if last_occurrence is not None:
        return lines[last_occurrence + 1 :]
    else:
        print("There are no lines found after the marker.")
        return None


# --- This method parses the JDFTx.out file --- #
def parse_output_file(file_path, verbose=False):
    energies = []
    oxidation_states = {}
    output_file_name = file_path.split("/")[-1].split(".")[0]
    electronic_scf = False
    vdw = False
    solvent = None
    file_path = (
        "Data-raw/" + file_path + ".out"
        if not file_path.endswith(".out")
        else "Data-raw/" + file_path
    )
    marker = "*************** JDFTx 1.7.0  ***************"
    lines = read_after_last_instance(file_path, marker)
    solvent_exists = True

    for line in lines:

        # Check for energy and iteration
        energy_match = re.search(
            r"IonicMinimize: Iter:\s+(\d+)\s+(?:F|Etot):\s+([-+]?\d*\.\d+|\d+)", line
        )
        if energy_match:
            iteration = int(energy_match.group(1))
            energy = float(energy_match.group(2))
            energies.append((iteration, energy))

        # Check for oxidation states
        if "# oxidation-state" in line:
            parts = line.split()
            element = parts[2]
            states = list(map(float, parts[3:]))
            oxidation_states[element] = states

        # Check for electronic scf
        if "scf" in line:
            electronic_scf = True

        # Check for solvent information
        if "fluid None" in line:
            solvent_exists = False
        elif "fluid-solvent" in line:
            solvent = line.split()[1]

        # Check for van der waals correction
        if 'van-der-waals' in line:
            vdw = True

    solvent = "Vacuum" if not solvent_exists else solvent
    # Create dataframes
    df_energies = pd.DataFrame(energies, columns=["Iteration", "Energy"])
    df_oxidation_states = pd.DataFrame(
        list(oxidation_states.items()), columns=["Element", "Oxidation States"]
    )

    # Print required information
    if verbose:
        print(df_energies)
        print(df_oxidation_states)
        print(f"Output file name: {output_file_name}")
        print(f"Electronic SCF used: {electronic_scf}")
        print(f"Solvent: {solvent}")

    if df_energies.empty:
        print(f"Error: df_energies is empty for {file_path}.")
        raise ValueError(f"Error: df_energies is empty for {file_path}.")
    if df_oxidation_states.empty:
        raise ValueError("Error: df_oxidation_states is empty.")

    return df_energies, df_oxidation_states, output_file_name, electronic_scf, solvent, vdw


# --- This method parses xsf files --- #
def parse_xsf_file(file_path: str, model: str, radius=6):
    file_path = (
        "Data-raw/" + file_path + ".xsf"
        if not file_path.endswith(".xsf")
        else "Data-raw/" + file_path
    )
    atomic_symbols = {
        1: "H",
        6: "C",
        7: "N",
        8: "O",
        11: "Na",
        16: "S",
        22: "Ti",
        26: "Fe",
        28: "Ni",
    }
    atomic_weights = {
        'H': 1.008,
        'C': 12.011,
        'N': 14.007,
        'O': 15.999,
        'Na': 22.990,
        'S': 32.066,
        'Ti': 47.867,
        'Fe': 55.845,
        'Ni': 58.693,
    }
    marker = "PRIMCOORD"
    lines = read_after_last_instance(file_path, marker)
    xyz = []
    for line in lines[1:]:
        parts = line.split()
        element = int(parts[0])
        x, y, z = map(float, parts[1:])
        xyz.append(
            {
                "element": atomic_symbols[element],
                "x": x,
                "y": y,
                "z": z,
            }
        )
    xyz_df = pd.DataFrame(xyz)
    # Separate into NaPS and Sorbent by extracting the NaPS
    NaPS_df = xyz_df[xyz_df["element"].isin(["Na", "S"])]
    # if NaPS is longer than 10 atoms (more than Na2S8), then need to filter
    if NaPS_df.shape[0] > 10:
        # get rows where the z column for the S element has a mode greater than 2
        S_rows = NaPS_df[NaPS_df["element"] == "S"]
        S_mode = S_rows["z"].mode()
        NaPS_df = NaPS_df[NaPS_df["z"] > max(S_mode)]
    if NaPS_df.shape[0] > 10:
        raise RuntimeError(f"NaPS has more than 10 atoms in {file_path}")
    # sorbent df is the negative of the NaPS df from the xyz df
    sorbent_df = xyz_df[~xyz_df.index.isin(NaPS_df.index)]
    # reset indices
    NaPS_df.reset_index(drop=True, inplace=True)
    sorbent_df.reset_index(drop=True, inplace=True)
    # Run x2t on the extracted NaPS and Sorbent in the adsorption file
    if not NaPS_df.empty:
        NaPS_bond_df, NaPS_angle_df = x2t.run_x2t(NaPS_df, args_filename=model, args_print=False, args_radius=radius)
        NaPS_bond_df["name"] = model
        NaPS_angle_df["name"] = model
        NaPS_bond_df = NaPS_bond_df[["name"] + [col for col in NaPS_bond_df.columns if col != "name"]]
        NaPS_angle_df = NaPS_angle_df[["name"] + [col for col in NaPS_angle_df.columns if col != "name"]]
    else:
        NaPS_bond_df = pd.DataFrame(columns=["name", "Atoms", "Mean", "Count", "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])
        NaPS_angle_df = pd.DataFrame(columns=["name", "Atoms", "Mean /°", "Count", "Median /°", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])
    if not sorbent_df.empty:
        sorbent_bond_df, sorbent_angle_df = x2t.run_x2t(sorbent_df, args_filename=model, args_print=False, args_radius=radius)
        sorbent_bond_df["name"] = model
        sorbent_angle_df["name"] = model
        sorbent_bond_df = sorbent_bond_df[["name"] + [col for col in sorbent_bond_df.columns if col != "name"]]
        sorbent_angle_df = sorbent_angle_df[["name"] + [col for col in sorbent_angle_df.columns if col != "name"]]
    else:
        sorbent_bond_df = pd.DataFrame(columns=["name", "Atoms", "Mean", "Count", "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])
        sorbent_angle_df = pd.DataFrame(columns=["name", "Atoms", "Mean /°", "Count", "Median /°", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])

    return NaPS_bond_df, NaPS_angle_df, sorbent_bond_df, sorbent_angle_df


# --- This method extracts the polysulfide from the string --- #
def extract_NaPS(filename: str):
    if any(substring in filename for substring in ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]):
        polysulfide = filename.split("@")[0]
        polysulfide = polysulfide.split("-")[0]
        return polysulfide
    else:
        return None


# --- This method extracts the sorbent from the string --- #
def extract_sorbent(filename: str):
    if any(substring in filename for substring in ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]) and "@" not in filename:
        return None
    else:
        sorbent = filename.split("-")[0]
        if "@" in sorbent:
            sorbent = sorbent.split("@")[1]
        if ".csv" in sorbent:
            sorbent = sorbent.split(".")[0]
        return sorbent


# --- This method compiles the energies and oxidation states of all sorbents in a subdirectory of Data-raw --- #
def extract_data(category: str, verbose=False):
    # get working directory
    working_directory = os.getcwd()

    # Get list of file names
    directory = working_directory + f"/Data-raw/{category}"
    file_list = os.listdir(directory)
    out_files = []
    xsf_files = []
    for file in file_list:
        if file.endswith(".out"):
            out_files.append(file)
        elif file.endswith(".xsf"):
            xsf_files.append(file)
        else:
            pass

    # Handle output files
    extracted_outputs = []
    for file in out_files:
        df_energies, df_oxidation_states, output_file_name, electronic_scf, solvent, vdw = (
            parse_output_file(f"{category}/{file}")
        )
        for _, row in df_energies.iterrows():
            extracted_outputs.append(
                {
                    "name": output_file_name,
                    "iteration": row["Iteration"],
                    "solvent": solvent,
                    "energy": row["Energy"],
                    "oxidation_states": df_oxidation_states.to_dict(orient="records"),
                    "electronic_scf": electronic_scf,
                    "vdw": vdw
                }
            )

    # create final output df        
    extracted_outputs_df = pd.DataFrame(extracted_outputs)
    extracted_outputs_df.sort_values(by="name", inplace=True)
    extracted_outputs_df["NaPS"] = extracted_outputs_df["name"].apply(extract_NaPS)
    extracted_outputs_df["sorbent"] = extracted_outputs_df["name"].apply(extract_sorbent)
    extracted_outputs_df.reset_index(drop=True, inplace=True)
    column_order = ["name", "NaPS", "sorbent", "solvent", "energy", "iteration", "electronic_scf", "vdw", "oxidation_states"]
    extracted_outputs_df = extracted_outputs_df[column_order]
    if verbose:
        print(f"df output for {category}")
        print(extracted_outputs_df.head())
    output_directory = os.path.join(working_directory, "Data-extracted/all")
    os.makedirs(output_directory, exist_ok=True)
    output_file_path = os.path.join(output_directory, f"{category}.csv")
    extracted_outputs_df.to_csv(output_file_path, index=False)
    final_extracted_outputs_df = extracted_outputs_df.loc[
        extracted_outputs_df.groupby("name")["iteration"].idxmax()
    ]
    output_directory = os.path.join(working_directory, "Data-extracted/final/energies")
    os.makedirs(output_directory, exist_ok=True)
    output_file_path = os.path.join(output_directory, f"{category}.csv")
    final_extracted_outputs_df.to_csv(output_file_path, index=False)

    # Handle xsf files
    NaPS_bond_dfs = []
    NaPS_angle_dfs = []
    sorbent_bond_dfs = []
    sorbent_angle_dfs = []

    for radius in [8, 62]:
        NaPS_extracted_xsf_bonds = []
        NaPS_extracted_xsf_angles = []
        sorbent_extracted_xsf_bonds = []
        sorbent_extracted_xsf_angles = []
        for file in xsf_files:
            NaPS_bond_df, NaPS_angle_df, sorbent_bond_df, sorbent_angle_df = parse_xsf_file(f"{category}/{file}", file.split(".")[0], radius=radius)
            # print(NaPS_bond_df)
            for _, row in NaPS_bond_df.iterrows():
                NaPS_extracted_xsf_bonds.append(row)
            for _, row in NaPS_angle_df.iterrows():
                NaPS_extracted_xsf_angles.append(row)
            for _, row in sorbent_bond_df.iterrows():
                sorbent_extracted_xsf_bonds.append(row)
            for _, row in sorbent_angle_df.iterrows():
                sorbent_extracted_xsf_angles.append(row)
        # Create Dataframes
        extracted_NaPS_xsf_bonds_df = pd.DataFrame(NaPS_extracted_xsf_bonds)
        extracted_NaPS_xsf_angles_df = pd.DataFrame(NaPS_extracted_xsf_angles)
        extracted_sorbent_xsf_bonds_df = pd.DataFrame(sorbent_extracted_xsf_bonds)
        extracted_sorbent_xsf_angles_df = pd.DataFrame(sorbent_extracted_xsf_angles)
        # print(extracted_NaPS_xsf_bonds_df)
        # Modify headers of dataframes
        if extracted_NaPS_xsf_bonds_df.empty:
            extracted_NaPS_xsf_bonds_df = pd.DataFrame(columns=["name", "NaPS", "sorbent", "Atoms", "Mean", "Count", "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])
            extracted_NaPS_xsf_angles_df = pd.DataFrame(columns=["name", "NaPS", "sorbent", "Atoms", "Mean", "Count", "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])
        if extracted_sorbent_xsf_bonds_df.empty:
            extracted_sorbent_xsf_bonds_df = pd.DataFrame(columns=["name", "NaPS", "sorbent", "Atoms", "Mean", "Count", "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])
            extracted_sorbent_xsf_angles_df = pd.DataFrame(columns=["name", "NaPS", "sorbent", "Atoms", "Mean", "Count", "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"])
        extracted_NaPS_xsf_bonds_df.sort_values(by="name", inplace=True)
        extracted_NaPS_xsf_angles_df.sort_values(by="name", inplace=True)
        extracted_sorbent_xsf_bonds_df.sort_values(by="name", inplace=True)
        extracted_sorbent_xsf_angles_df.sort_values(by="name", inplace=True)
        # get the NaPS column
        extracted_NaPS_xsf_bonds_df["NaPS"] = extracted_NaPS_xsf_bonds_df["name"].apply(extract_NaPS)
        extracted_NaPS_xsf_angles_df["NaPS"] = extracted_NaPS_xsf_angles_df["name"].apply(extract_NaPS)
        extracted_sorbent_xsf_bonds_df["sorbent"] = extracted_sorbent_xsf_bonds_df["name"].apply(extract_sorbent)
        extracted_sorbent_xsf_angles_df["sorbent"] = extracted_sorbent_xsf_angles_df["name"].apply(extract_sorbent)
        # get the sorbent column
        extracted_NaPS_xsf_bonds_df["sorbent"] = extracted_NaPS_xsf_bonds_df["name"].apply(extract_sorbent)
        extracted_NaPS_xsf_angles_df["sorbent"] = extracted_NaPS_xsf_angles_df["name"].apply(extract_sorbent)
        extracted_sorbent_xsf_bonds_df["NaPS"] = extracted_sorbent_xsf_bonds_df["name"].apply(extract_NaPS)
        extracted_sorbent_xsf_angles_df["NaPS"] = extracted_sorbent_xsf_angles_df["name"].apply(extract_NaPS)
        # reset index and combine with more columns from output df
        extracted_NaPS_xsf_bonds_df.reset_index(drop=True, inplace=True)
        extracted_NaPS_xsf_angles_df.reset_index(drop=True, inplace=True)
        extracted_sorbent_xsf_bonds_df.reset_index(drop=True, inplace=True)
        extracted_sorbent_xsf_angles_df.reset_index(drop=True, inplace=True)
        extracted_NaPS_xsf_bonds_df = pd.merge(extracted_NaPS_xsf_bonds_df, final_extracted_outputs_df[["name", "solvent", "electronic_scf", "vdw"]], on="name", how="left")
        extracted_NaPS_xsf_angles_df = pd.merge(extracted_NaPS_xsf_angles_df, final_extracted_outputs_df[["name", "solvent", "electronic_scf", "vdw"]], on="name", how="left")
        extracted_sorbent_xsf_bonds_df = pd.merge(extracted_sorbent_xsf_bonds_df, final_extracted_outputs_df[["name", "solvent", "electronic_scf", "vdw"]], on="name", how="left")
        extracted_sorbent_xsf_angles_df = pd.merge(extracted_sorbent_xsf_angles_df, final_extracted_outputs_df[["name", "solvent", "electronic_scf", "vdw"]], on="name", how="left")
        # reorder columns
        bond_column_order = ["name", "NaPS", "sorbent", "Atoms", "Mean", "solvent", "electronic_scf", "vdw", "Count", "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"]
        angles_column_order = ["name", "NaPS", "sorbent", "Atoms", "Mean", "solvent", "electronic_scf", "vdw", "Count",  "Median", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"]
        if not extracted_NaPS_xsf_bonds_df.empty:
            extracted_NaPS_xsf_bonds_df = extracted_NaPS_xsf_bonds_df[bond_column_order]
            extracted_NaPS_xsf_angles_df = extracted_NaPS_xsf_angles_df[angles_column_order]
        
        if not extracted_sorbent_xsf_bonds_df.empty:
            extracted_sorbent_xsf_bonds_df = extracted_sorbent_xsf_bonds_df[bond_column_order]
            extracted_sorbent_xsf_angles_df = extracted_sorbent_xsf_angles_df[angles_column_order]

        
        # append for final df creation
        NaPS_bond_dfs.append(extracted_NaPS_xsf_bonds_df)
        NaPS_angle_dfs.append(extracted_NaPS_xsf_angles_df)
        sorbent_bond_dfs.append(extracted_sorbent_xsf_bonds_df)
        sorbent_angle_dfs.append(extracted_sorbent_xsf_angles_df)

        bond_output_directory = os.path.join(working_directory, "Data-extracted/final/bonds")
        angle_output_directory = os.path.join(working_directory, "Data-extracted/final/angles")
        os.makedirs(bond_output_directory, exist_ok=True)
        os.makedirs(angle_output_directory, exist_ok=True)
        NaPS_bond_output_file_path = os.path.join(bond_output_directory, f"NaPS-{category}_bond-rad{radius}.csv")
        NaPS_angle_output_file_path = os.path.join(angle_output_directory, f"NaPS{category}_angle-rad{radius}.csv")
        sorbent_bond_output_file_path = os.path.join(bond_output_directory, f"sorbent-{category}_bond-rad{radius}.csv")
        sorbent_angle_output_file_path = os.path.join(angle_output_directory, f"sorbent-{category}_angle-rad{radius}.csv")
        extracted_NaPS_xsf_bonds_df.to_csv(NaPS_bond_output_file_path, index=False, encoding='utf-8-sig')
        extracted_NaPS_xsf_angles_df.to_csv(NaPS_angle_output_file_path, index=False, encoding='utf-8-sig')
        extracted_sorbent_xsf_bonds_df.to_csv(sorbent_bond_output_file_path, index=False, encoding='utf-8-sig')
        extracted_sorbent_xsf_angles_df.to_csv(sorbent_angle_output_file_path, index=False, encoding='utf-8-sig')

    print(f"Data extraction for {category} is complete.")
    return final_extracted_outputs_df, NaPS_bond_dfs, NaPS_angle_dfs, sorbent_bond_dfs, sorbent_angle_dfs


# --- This method extracts the adsorption energies using the finalized energies data --- #
def extract_adsorption_energies(sorbent: str, verbose=False):
    sorbent = sorbent if "NaPS@" not in sorbent else sorbent.split("@")[1]
    # get working directory
    working_directory = os.getcwd()
    directory = working_directory + f"/Data-extracted/final/energies"
    file_list = os.listdir(directory)
    # get NaPS & associated sorbent energies & combined model energies
    assert "NaPS.csv" in file_list, f"NaPS.csv not found in {directory}"
    assert f"Sorbents.csv" in file_list, f"Sorbents.csv not found in {directory}"
    assert f"NaPS@{sorbent}.csv" in file_list, f"NaPS@{sorbent}.csv not found in {directory}"
    # get dataframes
    NaPS_df = pd.read_csv(f"{directory}/NaPS.csv")
    sorbent_df = pd.read_csv(f"{directory}/Sorbents.csv")
    combined_df = pd.read_csv(f"{directory}/NaPS@{sorbent}.csv")
    
    # get binding energies
    adsorption_df = pd.DataFrame(["NaPS", "Sorbent", "Solvent", "Energy", "Electronic SCF", "vdw"])
    ads_energies = []
    for _, row in combined_df.iterrows():
        NaPS = row["NaPS"]
        sorbent = row["sorbent"]
        solvent = row["solvent"]
        electronic_scf = row["electronic_scf"]
        vdw = row["vdw"]
        # Get NaPS energy alone
        try:
            NaPS_energy = NaPS_df.loc[
                (NaPS_df["NaPS"] == NaPS) &
                (NaPS_df["solvent"] == solvent) &
                (NaPS_df["electronic_scf"] == electronic_scf) &
                (NaPS_df["vdw"] == vdw),
                "energy"
            ]
            if not NaPS_energy.empty:
                NaPS_energy = NaPS_energy.values[0]
            else:
                # NaPS_energy = None
                print(f"NaPS energy not found for {NaPS} in {solvent} in {sorbent} given electronic_scf: {electronic_scf} and vdw: {vdw}")
                NaPS_energy = 0
        except Exception as e:
            print(e)
        # Get sorbent energy alone
        try:
            sorbent_energy = sorbent_df.loc[
                (sorbent_df["sorbent"] == sorbent) &
                (sorbent_df["solvent"] == solvent) &
                (sorbent_df["electronic_scf"] == electronic_scf) &
                (sorbent_df["vdw"] == vdw),
                "energy"
            ]
            if not sorbent_energy.empty:
                sorbent_energy = sorbent_energy.values[0]
            else:
                # sorbent_energy = None
                print(f"Sorbent energy not found for {sorbent} in {solvent} given electronic_scf: {electronic_scf} and vdw: {vdw}")
                sorbent_energy = 0
        except Exception as e:
            print(e)
        # get combined energy
        sum_energy = sorbent_energy + NaPS_energy
        
        # get mixed model
        combined_energy = row["energy"]

        # adsorption energy
        ads_energy = sum_energy - combined_energy

        # combine rows to df
        ads_energies.append(
            {
                "NaPS": NaPS,
                "Sorbent": sorbent,
                "Solvent": solvent,
                "Binding Energy": ads_energy,
                "Combined Energy": combined_energy,
                "Sum Energy": sum_energy,
                "Electronic SCF": electronic_scf,
                "vdw": vdw
            }
        )

    adsorption_df = pd.DataFrame(ads_energies)
    output_directory = os.path.join(working_directory, "Data-extracted/final/binding_energies")
    os.makedirs(output_directory, exist_ok=True)
    output_file_path = os.path.join(output_directory, f"{sorbent}.csv")
    adsorption_df.to_csv(output_file_path, index=False)
    return adsorption_df

# Example usage
if __name__ == "__main__":
    # df_energies, df_oxidation_states, output_file_name, electronic_scf, solvent = (
    #     parse_output_file("Sorbents/FeN4-PC")
    # )
    # parse_xsf_file("Sorbents/FeN4-PC", "FeNe")
    # extract_data("Sorbents")
    # extract_data("NaPS@graphene_vdw")
    extract_data("NaPS@NiS2")
    # extract_adsorption_energies("NaPS@NiS2")
    # extract_adsorption_energies("NaPS@graphene_vdw")

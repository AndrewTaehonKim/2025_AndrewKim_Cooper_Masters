import pandas as pd
import os
import re

import xyz2tab as x2t


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
        raise ValueError("Error: df_energies is empty.")
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
    bond_df, angle_df = x2t.run_x2t(xyz_df, args_filename=model, args_print=False, args_radius=radius)
    bond_df["name"] = model
    angle_df["name"] = model
    bond_df = bond_df[["name"] + [col for col in bond_df.columns if col != "name"]]
    angle_df = angle_df[["name"] + [col for col in angle_df.columns if col != "name"]]
    return bond_df, angle_df


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
    bond_dfs = []
    angle_dfs = []
    for radius in [8, 62]:
        extracted_xsf_bonds = []
        extracted_xsf_angles = []
        for file in xsf_files:
            bond_df, angle_df = parse_xsf_file(f"{category}/{file}", file.split(".")[0], radius=radius)
            for _, row in bond_df.iterrows():
                extracted_xsf_bonds.append(row)
            for _, row in angle_df.iterrows():
                extracted_xsf_angles.append(row)
        
        extracted_xsf_bonds_df = pd.DataFrame(extracted_xsf_bonds)
        extracted_xsf_angles_df = pd.DataFrame(extracted_xsf_angles)
        extracted_xsf_bonds_df.sort_values(by="name", inplace=True)
        extracted_xsf_angles_df.sort_values(by="name", inplace=True)
        extracted_xsf_bonds_df["NaPS"] = extracted_xsf_bonds_df["name"].apply(extract_NaPS)
        extracted_xsf_angles_df["NaPS"] = extracted_xsf_angles_df["name"].apply(extract_NaPS)
        extracted_xsf_bonds_df["sorbent"] = extracted_xsf_bonds_df["name"].apply(extract_sorbent)
        extracted_xsf_angles_df["sorbent"] = extracted_xsf_angles_df["name"].apply(extract_sorbent)
        extracted_xsf_bonds_df.reset_index(drop=True, inplace=True)
        extracted_xsf_angles_df.reset_index(drop=True, inplace=True)
        extracted_xsf_bonds_df = extracted_xsf_bonds_df.merge(
            extracted_outputs_df[["name", "solvent", "electronic_scf", "vdw"]],
            on="name",
            how="left"
        )
        extracted_xsf_angles_df = extracted_xsf_angles_df.merge(
            extracted_outputs_df[["name", "solvent", "electronic_scf", "vdw"]],
            on="name",
            how="left"
        )
        bond_column_order = ["name", "NaPS", "sorbent", "Atoms", "solvent", "electronic_scf", "vdw", "Count", "Mean /Å", "Median /Å", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"]
        extracted_xsf_bonds_df = extracted_xsf_bonds_df[bond_column_order]
        angles_column_order = ["name", "NaPS", "sorbent", "Atoms", "solvent", "electronic_scf", "vdw", "Count", "Mean /°", "Median /°", "Sam. std. dev.", "Pop. std. dev.", "Std. error", "Skewness"]
        extracted_xsf_angles_df = extracted_xsf_angles_df[angles_column_order]
        bond_dfs.append(extracted_xsf_bonds_df)
        angle_dfs.append(extracted_xsf_angles_df)
        bond_output_directory = os.path.join(working_directory, "Data-extracted/final/bonds")
        angle_output_directory = os.path.join(working_directory, "Data-extracted/final/angles")
        os.makedirs(bond_output_directory, exist_ok=True)
        os.makedirs(angle_output_directory, exist_ok=True)
        bond_output_file_path = os.path.join(bond_output_directory, f"{category}_bond-rad{radius}.csv")
        angle_output_file_path = os.path.join(angle_output_directory, f"{category}_angle-rad{radius}.csv")
        extracted_xsf_bonds_df.to_csv(bond_output_file_path, index=False, encoding='utf-8-sig')
        extracted_xsf_angles_df.to_csv(angle_output_file_path, index=False, encoding='utf-8-sig')
    
    print(f"Data extraction for {category} is complete.")
    return final_extracted_outputs_df, bond_dfs[0], bond_dfs[1], angle_dfs[0], angle_dfs[1]


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
    # extract_data("NaPS@graphene")
    # extract_data("NaPS")
    extract_adsorption_energies("NaPS@NiS2")

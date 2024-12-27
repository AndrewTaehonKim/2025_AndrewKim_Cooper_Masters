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
        return []


# --- This method parses the JDFTx.out file --- #
def parse_output_file(file_path, print=False):
    energies = []
    oxidation_states = {}
    output_file_name = None
    electronic_scf = False
    solvent = None
    file_path = (
        "Data-raw/" + file_path + ".out"
        if not file_path.endswith(".out")
        else "Data-raw/" + file_path
    )
    marker = "*************** JDFTx 1.7.0  ***************"
    lines = read_after_last_instance(file_path, marker)

    for line in lines:
        # Check for energy and iteration
        energy_match = re.search(
            r"IonicMinimize: Iter:\s+(\d+)\s+F:\s+([-+]?\d*\.\d+|\d+)", line
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

        # Check for output file name
        if "dump-name" in line:
            output_file_name = line.split()[1].split(".")[0]

        # Check for electronic scf
        if "scf" in line:
            electronic_scf = True

        # Check for solvent information
        solvent_exists = True
        if "fluid None" in line:
            solvent_exists = False
        elif "fluid-solvent" in line:
            solvent = line.split()[1]

    # Create dataframes
    df_energies = pd.DataFrame(energies, columns=["Iteration", "Energy"])
    df_oxidation_states = pd.DataFrame(
        list(oxidation_states.items()), columns=["Element", "Oxidation States"]
    )

    # Print required information
    if print:
        print(df_energies)
        print(df_oxidation_states)
        print(f"Output file name: {output_file_name}")
        print(f"Electronic SCF used: {electronic_scf}")
        print(f"Solvent: {solvent if solvent_exists else 'None'}")

    return df_energies, df_oxidation_states, output_file_name, electronic_scf, solvent


# --- This method parses xsf files --- #
def parse_xsf_file(file_path):
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


# --- This method compiles the energies and oxidation states of all sorbents in a subdirectory of Data-raw --- #
def extract_data(category: str):
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
        df_energies, df_oxidation_states, output_file_name, electronic_scf, solvent = (
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
                }
            )
    extracted_outputs_df = pd.DataFrame(extracted_outputs)
    output_directory = os.path.join(working_directory, "Data-extracted")
    os.makedirs(output_directory, exist_ok=True)
    output_file_path = os.path.join(output_directory, f"{category}.csv")
    extracted_outputs_df.to_csv(output_file_path, index=False)

    # Handle xsf files


# Example usage
if __name__ == "__main__":
    # df_energies, df_oxidation_states, output_file_name, electronic_scf, solvent = (
    #     parse_output_file("Sorbents/FeN4-PC")
    # )
    parse_xsf_file("Sorbents/FeN4-PC")
    # extract_data("Sorbents")

import pandas as pd
import re

def parse_output_file(file_path):
    energies = []
    oxidation_states = {}
    output_file_name = None
    electronic_scf = False
    solvent = None
    file_path = 'Data-raw/' + file_path + '.out' if not file_path.endswith('.out') else 'Data-raw/' + file_path
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Check for energy and iteration
            energy_match = re.search(r'IonicMinimize: Iter:\s+(\d+)\s+F:\s+([-+]?\d*\.\d+|\d+)', line)
            if energy_match:
                iteration = int(energy_match.group(1))
                energy = float(energy_match.group(2))
                energies.append((iteration, energy))
            
            # Check for oxidation states
            if '# oxidation-state' in line:
                parts = line.split()
                element = parts[2]
                states = list(map(float, parts[3:]))
                oxidation_states[element] = states
            
            # Check for output file name
            if 'dump-name' in line:
                output_file_name = line.split()[1].split('.')[0]
            
            # Check for electronic scf
            if 'scf' in line:
                electronic_scf = True
            
            # Check for solvent information
            solvent_exists = True
            if 'fluid None' in line:
                solvent_exists = False
            elif 'fluid-solvent' in line:
                solvent = line.split()[1]

    # Create dataframes
    df_energies = pd.DataFrame(energies, columns=['Iteration', 'Energy'])
    df_oxidation_states = pd.DataFrame(list(oxidation_states.items()), columns=['Element', 'Oxidation States'])

    # Print required information
    print(f"Output file name: {output_file_name}")
    print(f"Electronic SCF used: {electronic_scf}")
    print(f"Solvent: {solvent if solvent_exists else 'None'}")

    return df_energies, df_oxidation_states

# Example usage
if __name__ == '__main__':
    df_energies, df_oxidation_states = parse_output_file('Sorbents/FeN4-PC')
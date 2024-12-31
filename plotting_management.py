import pandas as pd
import numpy as np
import os
import ast
import matplotlib.pyplot as plt
from data_management import extract_sorbent

import matplotlib
matplotlib.use('Agg')  # Use a non-GUI backend

dpi = 300

# --- This method plots the binding energies for each csv in Data-extracted/binding_energies --- #
def plot_binding_energies():
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/binding_energies')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv')]
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    solvents = ["Vacuum", "Glyme", "PC"]
    # Initialize a dictionary to hold the data
    # Read the data from the CSV files
    for csv in csvs:
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data = {solvent: [] for solvent in solvents}
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df = df[(df['NaPS'] == tick) & (df['Solvent'] == solvent)]
                if not filtered_df.empty:
                    data[solvent].append(filtered_df['Binding Energy'].values[0])
                else:
                    print(f"Missing energy data for {tick} in {solvent} in {csv}")
                    data[solvent].append(0)  # Handle missing data
        # Plot the data
        fig, ax = plt.subplots(figsize=(6, 5))
        bar_width = 0.2
        index = range(len(x_ticks))
        
        for i, solvent in enumerate(solvents):
            ax.bar([p + bar_width * i for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
        
        ax.set_xticks([p + bar_width for p in index])
        ax.set_xticklabels(x_ticks)
        y_min = min(min(values) for values in data.values())
        y_max = max(max(values) for values in data.values())
        y_range = y_max - y_min
        ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
        ax.set_ylabel('Binding Energy (hartrees)')
        ax.legend()
        
        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures/binding_energies')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        plt.savefig(os.path.join(figures_dir, f'{extract_sorbent(csv)}.jpg'), dpi=dpi)
        plt.close()


# --- This method plots the bond lengths of the NaPS molecules in the adsorbed state vs. non-adsorbed state --- #
def plot_NaPS_bonds():
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/bonds')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv')]
    NaPS_csvs = [csv for csv in csvs if '@' not in csv and "Sorbent" not in csv and "sorbent" not in csv]
    adsorption_csvs = [csv for csv in csvs if '@' in csv and 'NaPS-' in csv]
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    solvents = ["Vacuum", "Glyme", "PC"]
    # Extract unique rad# from the CSV filenames
    unique_rads = set()
    for csv in csvs:
        rad = csv.split('-')[-1]
        rad = rad.split('.')[0]
        unique_rads.add(rad)
    unique_rads = sorted(unique_rads)
    NaPS_rads = {rad: None for rad in unique_rads}
    # Get baseline bond lengths for NaPS
    for csv in NaPS_csvs:
        # save data for NaPS
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data_NaS = {solvent: [] for solvent in solvents}
        data_SS = {solvent: [] for solvent in solvents}
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df_NaS = df[((df['Atoms'] == 'S–Na')|(df['Atoms'] == 'Na–S')) & (df['NaPS'] == tick) & (df['solvent'] == solvent)]
                filtered_df_SS = df[(df['Atoms'] == 'S–S') & (df['NaPS'] == tick) & (df['solvent'] == solvent)]
                if not filtered_df_NaS.empty:
                    data_NaS[solvent].append(filtered_df_NaS['Mean'].values[0])
                else:
                    if solvent != 'PC':
                        print(f"Missing NaPS bond length data for {tick} in {solvent} in {csv} because Na-S is missing")
                    data_NaS[solvent].append(0)
                if not filtered_df_SS.empty:
                    data_SS[solvent].append(filtered_df_SS['Mean'].values[0])
                else:
                    if tick == "Na2S":
                        data_SS[solvent].append(0)
                    else:
                        if solvent != 'PC':
                            print(f"Missing NaPS bond length data for {tick} in {solvent} in {csv} because S-S is missing")
                        data_SS[solvent].append(0)
        # Add data to NaPS_rads based on the csv name
        rad = csv.split('-')[-1].split('.')[0]
        NaPS_rads[rad] = {'NaS': data_NaS, 'SS': data_SS}
    # Read the data from the CSV files
    for csv in adsorption_csvs:
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data_NaS = {solvent: [] for solvent in solvents}
        data_SS = {solvent: [] for solvent in solvents}
        rad = csv.split('-')[-1].split('.')[0]
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df_NaS = df[((df['Atoms'] == 'S–Na')|(df['Atoms'] == 'Na–S')) & (df['NaPS'] == tick) & (df['solvent'] == solvent)]
                filtered_df_SS = df[(df['Atoms'] == 'S–S') & (df['NaPS'] == tick) & (df['solvent'] == solvent)]
                if not filtered_df_NaS.empty:
                    data_NaS[solvent].append(filtered_df_NaS['Mean'].values[0])
                else:
                    if solvent != 'PC':
                        print(f"Missing NaPS bond length data for {tick} in {solvent} in {csv} because Na-S is missing")
                    data_NaS[solvent].append(0)
                if not filtered_df_SS.empty:
                    data_SS[solvent].append(filtered_df_SS['Mean'].values[0])
                else:
                    if tick == "Na2S":
                        data_SS[solvent].append(0)
                    else:
                        if solvent != 'PC':
                            print(f"Missing NaPS bond length data for {tick} in {solvent} in {csv} because S-S is missing")
                        data_SS[solvent].append(0)
        # Plot the data as two subplots in one figure, the left for Na-S bonds and the right for S-S bonds
        fig, (ax_NaS, ax_SS) = plt.subplots(1, 2, figsize=(12, 10))
        bar_width = 0.2
        NaS_index = range(len(x_ticks))
        SS_index = range(len(x_ticks[1:]))

        # Plot bars
        for i, solvent in enumerate(solvents):
            ax_NaS.bar([p + bar_width * i for p in NaS_index], data_NaS[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
            ax_SS.bar([p + bar_width * i for p in SS_index], data_SS[solvent][1:], bar_width, label=solvent, hatch=['', '\\', '/'][i])
            # Plot thick horizontal line for baseline NaS and SS bond lengths using NaPS_rads
            base_NaS = NaPS_rads[rad]['NaS']
            base_SS = NaPS_rads[rad]['SS']
            for j in NaS_index:
                if base_NaS[solvent][j] != 0:
                    if j == 0 and solvent == 'Vacuum':
                        ax_NaS.hlines(base_NaS[solvent][j], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3, label='Isolated NaPS')
                    else:
                        ax_NaS.hlines(base_NaS[solvent][j], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3)
            for j in SS_index:
                if base_SS[solvent][j+1] != 0:
                    if j == 0 and solvent == 'Vacuum':
                        ax_SS.hlines(base_SS[solvent][j+1], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3, label='Isolated NaPS')
                    else:
                        ax_SS.hlines(base_SS[solvent][j+1], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3)
        
        
        # Set x-ticks and x-tick labels
        ax_NaS.set_xticks([p + bar_width for p in NaS_index])
        ax_NaS.set_xticklabels(x_ticks)
        ax_SS.set_xticks([p + bar_width for p in SS_index])
        ax_SS.set_xticklabels(x_ticks[1:])

        # Calculate y-axis limits
        y_min_NaS = min(min(values) for values in data_NaS.values())
        y_max_NaS = max(max(values) for values in data_NaS.values())
        y_range_NaS = y_max_NaS - y_min_NaS
        ax_NaS.set_ylim(y_min_NaS - 0.1 * y_range_NaS if y_min_NaS > 0 else 1.5, y_max_NaS + 0.1 * y_range_NaS)

        y_min_SS = min(min(values) for values in data_SS.values())
        y_max_SS = max(max(values) for values in data_SS.values())
        y_range_SS = y_max_SS - y_min_SS
        ax_SS.set_ylim(y_min_SS - 0.1 * y_range_SS  if y_min_SS > 0 else 1.5, y_max_SS + 0.1 * y_range_SS)

        # Set y-axis labels and legends
        ax_NaS.set_ylabel('Bond Length (Å)')
        ax_SS.set_ylabel('Bond Length (Å)')
        ax_NaS.legend()
        ax_SS.legend()

        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures/NaPS-bond_lengths')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        plt.savefig(os.path.join(figures_dir, f'{csv.split(".")[0]}.jpg'), bbox_inches='tight', pad_inches=0.1, dpi=dpi)
        plt.close()

# --- This method plots the adsorption bond lengths between the NaPS molecules and the sorbent --- #
def plot_adsorption_lengths():
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/adsorption_lengths')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv')]
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    solvents = ["Vacuum", "Glyme", "PC"]
    # Initialize a dictionary to hold the data
    # Read the data from the CSV files
    for csv in csvs:
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data = {solvent: [] for solvent in solvents}
        connections = {solvent: [] for solvent in solvents}
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df = df[(df['NaPS'] == tick) & (df['solvent'] == solvent)]
                if not filtered_df.empty:
                    data[solvent].append(filtered_df['Distance'].values[0])
                    connections[solvent].append(
                        str(filtered_df['NaPS Atom'].values[0]) + '-' + str(filtered_df['Sorbent Atom'].values[0])
                    )
                else:
                    print(f"Missing adsorption bond length data for {tick} in {solvent} in {csv}")
                    data[solvent].append(0)
                    connections[solvent].append('')
        # Plot the data
        fig, ax = plt.subplots(figsize=(6, 5))
        bar_width = 0.2
        index = range(len(x_ticks))

        for i, solvent in enumerate(solvents):
            ax.bar([p + bar_width * i for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
            # label each bar with the NaPS Atom - Sorbent Atom
            for j, value in enumerate(data[solvent]):
                ax.text(j + bar_width * i, value + 0.05, connections[solvent][j], ha='center', va='bottom', size=5)

        ax.set_xticks([p + bar_width for p in index])
        ax.set_xticklabels(x_ticks)
        y_min = min(min(values) for values in data.values())
        y_max = max(max(values) for values in data.values())
        y_range = y_max - y_min
        ax.set_ylim(1.5, y_max + 0.1 * y_range)
        ax.set_ylabel('Adsorption Distance (Å)')
        ax.legend()
        
        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures/adsorption_lengths')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        plt.savefig(os.path.join(figures_dir, f'{extract_sorbent(csv)}.jpg'), dpi=dpi)
        plt.close()


# --- This method plots the average oxidation states for Na and S in the NaPS molecules --- #
def plot_oxidation_states():
# Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/energies')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv')]
    NaPS_csvs = [csv for csv in csvs if '@' not in csv and "Sorbent" not in csv and "sorbent" not in csv]
    adsorption_csvs = [csv for csv in csvs if '@' in csv]
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    solvents = ["Vacuum", "Glyme", "PC"]
    # Extract unique rad# from the CSV filenames
    # Get baseline bond lengths for NaPS
    NaPS_data = {}
    for csv in NaPS_csvs:
        # save data for NaPS
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data_Na = {solvent: [] for solvent in solvents}
        data_S = {solvent: [] for solvent in solvents}
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df = df[(df['NaPS'] == tick) & (df['solvent'] == solvent)]
                if not filtered_df.empty:
                    data_Na[solvent].append(np.mean(np.array(ast.literal_eval(filtered_df['oxidation_states'].values[0])['Na'])))
                    data_S[solvent].append(np.mean(np.array(ast.literal_eval(filtered_df['oxidation_states'].values[0])['S'])))
                else:
                    data_Na[solvent].append(0)
                    data_S[solvent].append(0)
    NaPS_data = {'Na': data_Na, 'S': data_S}
    # Read the data from the CSV files
    for csv in adsorption_csvs:
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data_Na = {solvent: [] for solvent in solvents}
        data_S = {solvent: [] for solvent in solvents}
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df = df[(df['NaPS'] == tick) & (df['solvent'] == solvent)]
                if not filtered_df.empty:
                    data_Na[solvent].append(np.mean(np.array(ast.literal_eval(filtered_df['oxidation_states'].values[0])['Na'])))
                    data_S[solvent].append(np.mean(np.array(ast.literal_eval(filtered_df['oxidation_states'].values[0])['S'])))
                else:
                    data_Na[solvent].append(0)
                    data_S[solvent].append(0)
        # Plot the data as two subplots in one figure, the left for Na-S bonds and the right for S-S bonds
        fig, (ax_Na, ax_S) = plt.subplots(1, 2, figsize=(12, 10))
        bar_width = 0.2
        index = range(len(x_ticks))
        # Plot bars
        for i, solvent in enumerate(solvents):
            ax_Na.bar([p + bar_width * i for p in index], data_Na[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
            ax_S.bar([p + bar_width * i for p in index], data_S[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
            # Plot thick horizontal line for baseline NaS and SS bond lengths using NaPS_rads
            base_Na = NaPS_data['Na'][solvent]
            base_S = NaPS_data['S'][solvent]
            for j in index:
                if base_Na[j] != 0:
                    if j == 0 and solvent == 'Vacuum':
                        ax_Na.hlines(base_Na[j], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3, label='Isolated NaPS')
                    else:
                        ax_Na.hlines(base_Na[j], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3)
                if base_S[j] != 0:
                    if j == 0 and solvent == 'Vacuum':
                        ax_S.hlines(base_S[j], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3, label='Isolated NaPS')
                    else:
                        ax_S.hlines(base_S[j], j + bar_width * i - bar_width / 2, j + bar_width * i + bar_width / 2, colors='red', linewidth=4, zorder=3)
        
        # Set x-ticks and x-tick labels
        ax_Na.set_xticks([p + bar_width for p in index])
        ax_Na.set_xticklabels(x_ticks)
        ax_S.set_xticks([p + bar_width for p in index])
        ax_S.set_xticklabels(x_ticks)

        # Calculate y-axis limits
        # y_min_Na = min(min(values) for values in data_Na.values())
        # y_max_Na = max(max(values) for values in data_Na.values())
        # y_range_Na = y_max_Na - y_min_Na
        # ax_Na.set_ylim(y_min_Na - 0.1 * y_range_Na, y_max_Na + 0.1 * y_range_Na)

        # y_min_S = min(min(values) for values in data_S.values())
        # y_max_S = max(max(values) for values in data_S.values())
        # y_range_S = y_max_S - y_min_S
        # ax_S.set_ylim(y_min_S - 0.1 * y_range_S, y_max_S + 0.1 * y_range_S)

        # Set y-axis labels and legends
        ax_Na.set_ylabel('Bond Length (Å)')
        ax_S.set_ylabel('Bond Length (Å)')
        ax_Na.legend()
        ax_S.legend()

        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures/oxidation_states')
        os.makedirs(figures_dir, exist_ok=True)
        plt.savefig(os.path.join(figures_dir, f'{csv.split(".")[0]}.jpg'), bbox_inches='tight', pad_inches=0.1, dpi=dpi)


if __name__ == "__main__":
    # plot_binding_energies()
    # plot_NaPS_bonds()
    # plot_adsorption_lengths()
    plot_oxidation_states()
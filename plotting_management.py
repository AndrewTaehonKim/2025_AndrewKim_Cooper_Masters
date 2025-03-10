import pandas as pd
import numpy as np
import os
import ast
import matplotlib.pyplot as plt
from data_management import extract_sorbent
from matplotlib.ticker import FuncFormatter
import matplotlib
matplotlib.use('Agg')  # Use a non-GUI backend
matplotlib.rcParams['text.usetex'] = False

dpi = 300
conversion_factor = 2625.5  # Convert hartree to kJ/mol
NaPS_labels = ["$Na_2S$", "$Na_2S_2$", "$Na_2S_4$", "$Na_2S_6$", "$Na_2S_8$"]

# --- This method plots the convergence of NaPSs based on cell size --- #
def plot_cell_size_convergence():
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/convergence')
    csv = "NaPS_cell_convergence"
    csv_path = os.path.join(directory, csv + ".csv")
    df = pd.read_csv(csv_path)
    # get unique cell sizes
    cell_sizes = df['Cell Length'].unique()
    NaPSs = df['NaPS'].unique()
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    for i, NaPS in enumerate(NaPSs):
        filtered_df = df[df['NaPS'] == NaPS]
        filtered_df.dropna(inplace=True)
        filtered_df.drop(filtered_df[filtered_df['%Difference'] == 0].index, inplace=True)
        indices = [i for i in range(filtered_df.shape[0])]
        axs[i//3][i%3].set_title(f'{NaPS_labels[i]}')
        axs[i//3][i%3].scatter(indices, np.abs(filtered_df['%Difference'].to_numpy()))
        axs[i//3][i%3].set_xticks(indices)
        axs[i//3][i%3].set_xticklabels(np.array(filtered_df['Cell Length']))#, rotation=45)
        axs[i//3][i%3].set_xlabel('Cell Length (Å)')
        axs[i//3][i%3].set_ylabel('Absolute % Difference (%)')
        axs[i//3][i%3].set_ylim(0, 1.1e-5)
    plt.subplots_adjust(wspace=.4, hspace=.7)
    plt.savefig(os.path.join(os.getcwd(), 'Figures/NaPS-cell_size_convergence.jpg'), dpi=dpi*2, bbox_inches='tight')

# --- This method plots the kpoint convergence and time taken --- #
def plot_kpoint_convergence():
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/convergence')
    csv = "graphene_kpoint_convergence"
    csv_path = os.path.join(directory, csv + ".csv")
    df = pd.read_csv(csv_path)
    # get unique kpoints
    kpoints = df['kpoint'].unique()
    fig, axs = plt.subplots(1, 2, figsize=(8, 2))
    # plot kpoint energy convergence
    axs[0].bar(kpoints, df['energy'])
    axs[0].set_xlabel('kpoint folding')
    axs[0].set_ylabel('Energy (Hartree)')
    axs[0].set_ylim(-417, -410)
    # plot time required
    axs[1].bar(kpoints, df['Time (min)'])
    axs[1].set_xlabel('kpoint folding')
    axs[1].set_ylabel('Time (min)')    
    
    plt.savefig(os.path.join(os.getcwd(), 'Figures/graphene_kpoint_convergence.jpg'), dpi=dpi, bbox_inches='tight')


# --- This method plots the energyies for the NaPSs with vdw, without vdw, scf, and without scf --- #
def plot_NaPS_energy():
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/energies')
    csv = "NaPS.csv"
    csv_path = os.path.join(directory, csv)
    full_df = pd.read_csv(csv_path)
    # create a subplot of 2 by 2 with SCF and elecmin as first and second row and no vdw and vdw as first and second column
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    x_labels = ["$Na_2S$", "$Na_2S_2$", "$Na_2S_4$", "$Na_2S_6$", "$Na_2S_8$"]
    solvents = ["Vacuum", "Glyme", "PC"]
    # Initialize a dictionary to hold the data
    nvdw_scf_data = {solvent: [] for solvent in solvents}
    vdw_scf_data = {solvent: [] for solvent in solvents}
    nvdw_elecmin_data = {solvent: [] for solvent in solvents}
    vdw_elecmin_data = {solvent: [] for solvent in solvents}

    # split data into 4
    nvdw_scf_df = full_df[(full_df['vdw'] == False) & (full_df['electronic_scf'] == True)]
    vdw_scf_df = full_df[(full_df['vdw'] == True) & (full_df['electronic_scf'] == True)]
    nvdw_elecmin_df = full_df[(full_df['vdw'] == False) & (full_df['electronic_scf'] == False)]
    vdw_elecmin_df = full_df[(full_df['vdw'] == True) & (full_df['electronic_scf'] == False)]

    # Read the data from the CSV files
    for df, data in zip([nvdw_scf_df, vdw_scf_df, nvdw_elecmin_df, vdw_elecmin_df], [nvdw_scf_data, vdw_scf_data, nvdw_elecmin_data, vdw_elecmin_data]):
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df = df[(df['NaPS'] == tick) & (df['solvent'] == solvent)]
                if not filtered_df.empty:
                    data[solvent].append(filtered_df['energy'].values[0]*-1)
                else:
                    print(f"Missing energy data for {tick} in {solvent} in {csv}")
                    data[solvent].append(0)
    # Plot the data
    bar_width = 0.2
    index = range(len(x_ticks))
    for i, data in enumerate([nvdw_scf_data, vdw_scf_data, nvdw_elecmin_data, vdw_elecmin_data]):
        for j, solvent in enumerate(solvents):
            axs[i // 2, i % 2].bar([p + bar_width * j for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][j])
        axs[i // 2, i % 2].set_xticks([p + bar_width for p in index])
        axs[i // 2, i % 2].set_xticklabels(x_ticks)
        axs[i // 2, i % 2].set_ylabel('Energy (Hartree)')
        axs[i // 2, i % 2].set_ylim(bottom=100)
        axs[i // 2, i % 2].legend()
    plt.savefig(os.path.join(os.getcwd(), 'Figures/NaPS-4_subplot-system_type.jpg'), dpi=dpi, bbox_inches='tight')

    # split data by NaPS
    Na2S_df = full_df[full_df['NaPS'] == 'Na2S']
    Na2S2_df = full_df[full_df['NaPS'] == 'Na2S2']
    Na2S4_df = full_df[full_df['NaPS'] == 'Na2S4']
    Na2S6_df = full_df[full_df['NaPS'] == 'Na2S6']
    Na2S8_df = full_df[full_df['NaPS'] == 'Na2S8']
    NaPS_df = [Na2S_df, Na2S2_df, Na2S4_df, Na2S6_df, Na2S8_df]

    # dictionaries to store
    Na2S_data = {solvent: [] for solvent in solvents}
    Na2S2_data = {solvent: [] for solvent in solvents}
    Na2S4_data = {solvent: [] for solvent in solvents}
    Na2S6_data = {solvent: [] for solvent in solvents}
    Na2S8_data = {solvent: [] for solvent in solvents}
    NaPS_data = [Na2S_data, Na2S2_data, Na2S4_data, Na2S6_data, Na2S8_data]
    # plotting variables
    NaPSs = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    x_ticks = ["no vdw / scf", "vdw / scf", "no vdw / elecmin", "vdw / elecmin"]
    x_ticks = ["A", "B", "C", "D"]
    system_types = [(False, True), (True, True), (False, False), (True, False)] # vdw, scf
    # Plot the data
    fig, axs = plt.subplots(5, 1, figsize=(10, 15))
    for df, data in zip(NaPS_df, NaPS_data):
        for solvent in solvents:
            for i, tick in enumerate(x_ticks):
                filtered_df = df[(df['solvent'] == solvent) & (df['electronic_scf'] == system_types[i][0]) & (df['vdw'] == system_types[i][1])]
                if not filtered_df.empty:
                    data[solvent].append(filtered_df['energy'].values[0]*-1)
                else:
                    print(f"Missing energy data for {tick} in {solvent} in {csv}")
                    data[solvent].append(0)
    index = range(len(x_ticks))
    bar_width = 0.2
    for i, data in enumerate(NaPS_data):
        y_min = min(min(values) for values in data.values())
        y_max = max(max(values) for values in data.values())
        for j, solvent in enumerate(solvents):
            axs[i].bar([p + bar_width * j for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][j])
        axs[i].set_xticks([p + bar_width for p in index])
        axs[i].set_xticklabels(x_ticks)
        axs[i].set_ylabel('Energy (Hartree)')
        axs[i].set_title(f'{x_labels[i]}')
        axs[i].set_ylim(y_min-0.1, y_max+0.1)
        axs[i].legend()
    plt.savefig(os.path.join(os.getcwd(), 'Figures/NaPS-5_subplot-NaPS.jpg'), dpi=dpi, bbox_inches='tight')


# --- This method plots the oxidation states of the Na and S in the NaPS molecules --- #
def plot_NaPS_oxidation():
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/energies')
    csv = "NaPS.csv"
    csv_path = os.path.join(directory, csv)
    full_df = pd.read_csv(csv_path)
    # create a subplot of 2 by 2 with SCF and elecmin as first and second row and no vdw and vdw as first and second column
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    x_labels = ["$Na_2S$", "$Na_2S_2$", "$Na_2S_4$", "$Na_2S_6$", "$Na_2S_8$"]
    solvents = ["Vacuum", "Glyme", "PC"]
    # Initialize a dictionary to hold the data
    data_S = {solvent: [] for solvent in solvents}
    data_Na = {solvent: [] for solvent in solvents}
    data_S_std = {solvent: [] for solvent in solvents}
    data_Na_std = {solvent: [] for solvent in solvents}
    # Read the data from the CSV files
    for solvent in solvents:
        for tick in x_ticks:
            oxidation_list_S = []
            oxidation_list_Na = []
            filtered_df = full_df[(full_df['NaPS'] == tick) & (full_df['solvent'] == solvent)]
            if not filtered_df.empty:
                oxidation_list_Na.append(ast.literal_eval(filtered_df['oxidation_states'].values[0])['Na'])
                oxidation_list_S.append(ast.literal_eval(filtered_df['oxidation_states'].values[0])['S'])
            else:
                print(f"Missing oxidation state data for {tick} in {solvent} in {csv}")
                oxidation_list_Na.append(0)
                oxidation_list_S.append(0)
            data_Na[solvent].append(np.mean(np.array(oxidation_list_Na)))
            data_S[solvent].append(np.mean(np.array(oxidation_list_S)))
            data_Na_std[solvent].append(np.std(np.array(oxidation_list_Na)))
            data_S_std[solvent].append(np.std(np.array(oxidation_list_S)))
    # Plot the data
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    bar_width = 0.2
    index = range(len(x_ticks))
    for i, solvent in enumerate(solvents):
        # bar plot with std
        axs[0].bar([p + bar_width * i for p in index], data_Na[solvent], bar_width, label=solvent, yerr=data_Na_std[solvent], hatch=['', '\\', '/'][i], error_kw={'ecolor': 'black', 'capsize': 5, 'elinewidth': 1, 'capthick': 1})
        axs[1].bar([p + bar_width * i for p in index], data_S[solvent], bar_width, label=solvent, yerr=data_S_std[solvent], hatch=['', '\\', '/'][i], error_kw={'ecolor': 'black', 'capsize': 5, 'elinewidth': 1, 'capthick': 1})
    axs[0].set_xticks([p + bar_width for p in index])
    axs[0].set_xticklabels(x_labels)
    axs[1].set_xticks([p + bar_width for p in index])
    axs[1].set_xticklabels(x_labels)
    axs[0].set_ylabel('Oxidation State')
    axs[1].set_ylabel('Oxidation State')
    axs[0].legend()
    axs[1].legend()
    axs[0].set_ylim(0.6, 1.1)
    
    axs[0].set_title('Na Oxidation States')
    axs[1].set_title('S Oxidation States')
    plt.savefig(os.path.join(os.getcwd(), 'Figures/oxidation_states/NaPS.jpg'), dpi=dpi, bbox_inches='tight')

# --- This method plots the binding energies for each csv in Data-extracted/binding_energies --- #
"""
def plot_binding_energies():
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/binding_energies')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv') and 'reference' not in csv]
    
    # Get reference data
    reference_df = pd.read_csv(os.path.join(directory, 'references.csv'))
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    
    # correct for TiO2

    solvents = ["Vacuum", "Glyme", "PC"]
    # Initialize a dictionary to hold the data
    # Read the data from the CSV files
    for csv in csvs:
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data = {solvent: [] for solvent in solvents}
        sorbent = df['Sorbent'].values[0]
        # Save absolute minimum and maximum for plotting purposes
        y_min = 0
        y_max = 0
        for solvent in solvents:
            for tick in x_ticks:
                filtered_df = df[(df['NaPS'] == tick) & (df['Solvent'] == solvent)]
                if not filtered_df.empty:
                    # convert hartree to kJ/mol
                    data[solvent].append(filtered_df['Binding Energy'].values[0]*conversion_factor)
                else:
                    print(f"Missing energy data for {tick} in {solvent} in {csv}")
                    data[solvent].append(0)  # Handle missing data
        y_min = min(min(values) for values in data.values())
        y_max = max(max(values) for values in data.values())                    
        # Plot the data
        fig, ax = plt.subplots(figsize=(6, 5))
        bar_width = 0.2
        index = range(len(x_ticks))
        
        for i, solvent in enumerate(solvents):
            ax.bar([p + bar_width * i for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
        
        # Plot thick horizontal line for reference data based on matching NaPS and sorbent
        labeled = False
        energy_list = []
        energy_std_list = []
        for j in index:
            reference_data = reference_df[(reference_df['NaPS'] == x_ticks[j]) & (reference_df['Sorbent'] == sorbent)]
            # take average and std
            energy = reference_data['Binding Energy (kJ/mol)'].mean()
            energy_std = reference_data['Binding Energy (kJ/mol)'].std()
            energy_min = reference_data['Binding Energy (kJ/mol)'].min()
            energy_max = reference_data['Binding Energy (kJ/mol)'].max()
            if energy_min < y_min:
                y_min = energy_min
            if energy_max > y_max:
                y_max = energy_max
            # Convert the average, minimum, and maximum binding energy to kJ/mol
            energy_min_kJ = energy_min
            energy_max_kJ = energy_max
            # Calculate the asymmetric error bars
            yerr = [[energy - energy_min_kJ], [energy_max_kJ - energy]]
            # if reference data is larger than 1, print reference data, energy, and std
            # if reference_data.shape[0] > 1:
            #     print(f"Reference data for {sorbent} and {x_ticks[j]}: {reference_data['Binding Energy (kJ/mol)'].values}")
            #     print(f"Energy: {energy}")
            #     print(f"Std: {energy_std}")
            if not reference_data.empty:
                if labeled is False:
                    # label with a horizontal line
                    ax.hlines(energy, j - bar_width / 2, j + bar_width / 2, colors='red', linewidth=4, zorder=3, label=f'Reference')
                    labeled = True
                else:
                    ax.hlines(energy, j - bar_width / 2, j + bar_width / 2, colors='red', linewidth=4, zorder=3)
                if not pd.isna(energy_std):
                    # ax.errorbar(j, energy, yerr=energy_std, fmt='o', color='red', zorder=3, capsize=5, capthick=1, elinewidth=1)
                    ax.errorbar(j, energy, yerr=yerr, fmt='o', color='red', zorder=3, capsize=5, capthick=1, elinewidth=1)
        ax.set_xticks([p + bar_width for p in index])
        ax.set_xticklabels(NaPS_labels)
        # y_min = min(min(values) for values in data.values())
        # y_max = max(max(values) for values in data.values())
        print(f"y_min: {y_min}, y_max: {y_max}")
        y_max = y_max if y_max < 5000 else 500
        y_range = y_max - y_min
        ax.set_ylim((y_min - 0.1 * y_range) if y_min < 0 else 0, y_max + 0.1 * y_range)
        ax.set_ylabel('Binding Energy ($kJ\;mol^{-1}$)')
        ax.legend()
        
        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures/binding_energies')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        plt.savefig(os.path.join(figures_dir, f'{extract_sorbent(csv)}.jpg'), dpi=dpi)
        plt.close()
"""
def plot_binding_energies():
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/binding_energies')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv') and 'reference' not in csv]
    
    # Get reference data
    reference_df = pd.read_csv(os.path.join(directory, 'references.csv'))
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    solvents = ["Reference", "Vacuum", "Glyme", "PC"]
    # Initialize a dictionary to hold the data
    # Read the data from the CSV files
    for csv in csvs:
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        data = {solvent: [] for solvent in solvents}
        error_bars = {solvent: [] for solvent in solvents}
        sorbent = df['Sorbent'].values[0]
        for solvent in solvents[1:]:  # Skip "Reference" for now
            for tick in x_ticks:
                filtered_df = df[(df['NaPS'] == tick) & (df['Solvent'] == solvent)]
                if not filtered_df.empty:
                    # convert hartree to kJ/mol
                    data[solvent].append(filtered_df['Binding Energy'].values[0] * conversion_factor)
                else:
                    print(f"Missing energy data for {tick} in {solvent} in {csv}")
                    data[solvent].append(0)  # Handle missing data
                    error_bars[solvent].append([[0], [0]])  # No error bars
        # Add reference data
        reference_exists = True
        for tick in x_ticks:
            reference_data = reference_df[(reference_df['NaPS'] == tick) & (reference_df['Sorbent'] == sorbent)]
            if not reference_data.empty:
                # Calculate the average, minimum, and maximum binding energy
                average_binding_energy_kJ = reference_data['Binding Energy (kJ/mol)'].mean()
                energy_min_kJ = reference_data['Binding Energy (kJ/mol)'].min()
                energy_max_kJ = reference_data['Binding Energy (kJ/mol)'].max()
                data["Reference"].append(average_binding_energy_kJ)
                # Store the error bars
                yerr = [[average_binding_energy_kJ - energy_min_kJ], [energy_max_kJ - average_binding_energy_kJ]]
                error_bars["Reference"].append(yerr)
                print(f"Reference data for {sorbent} and {tick}: {average_binding_energy_kJ}")
            else:
                print(f"Missing reference data for {tick} in {sorbent}")
                reference_exists = False

        # Plot the data
        fig, ax = plt.subplots(figsize=(10, 6))
        bar_width = 0.2
        index = range(len(x_ticks))
        
        if reference_exists:
            for i, solvent in enumerate(solvents):
                if solvent == "Reference":
                    print(f"For sorbent {sorbent} ; Reference data: {data[solvent]}")
                    for j, tick in enumerate(x_ticks):
                        ax.bar(j + bar_width * i, data[solvent][j], bar_width, label="Literature" if j == 0 else "", color='red', yerr=error_bars[solvent][j], error_kw={'ecolor': 'black', 'capsize': 5, 'elinewidth': 2, 'capthick': 2})
                else:
                    ax.bar([p + bar_width * i for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i-1])
        else:
            for i, solvent in enumerate(solvents[1:]):  # Skip "Reference"
                ax.bar([p + bar_width * i for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
        
        ax.set_xticks([p + bar_width * (1.5 if reference_exists else 1) for p in index])
        x_tick_labels = ["$Na_2S$", "$Na_2S_2$", "$Na_2S_4$", "$Na_2S_6$", "$Na_2S_8$"]
        ax.set_xticklabels(x_tick_labels)
        
        # Filter out empty lists before calculating y_min and y_max
        non_empty_data = [values for values in data.values() if values]
        y_min = min(min(values) for values in non_empty_data)
        y_max = max(max(values) for values in non_empty_data)
        # if reference exists, check for y_max in error bars
        if reference_exists:
            y_max = max(y_max, np.max([error_bar[1][0]+data["Reference"] for error_bar in error_bars["Reference"]]))
        y_range = y_max - y_min
        ax.set_ylim((y_min - 0.1 * y_range) if y_min < 0 else 0, y_max + 0.1 * y_range)
        ax.set_ylabel('Binding Energy ($kJ\\;mol^{-1}$)', fontsize=15)

        ax.set_xlabel('NaPS', fontsize=15)
        ax.legend()
        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures/binding_energies')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        plt.savefig(os.path.join(figures_dir, f'{sorbent}.jpg'), dpi=300)
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
    x_tick_labels = [r"$Na_2S$", r"$Na_2S_2$", r"$Na_2S_4$", r"$Na_2S_6$", r"$Na_2S_8$"]
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
        fig, ax = plt.subplots(figsize=(10, 6))
        bar_width = 0.2
        index = range(len(x_ticks))

        for i, solvent in enumerate(solvents):
            ax.bar([p + bar_width * i for p in index], data[solvent], bar_width, label=solvent, hatch=['', '\\', '/'][i])
            # label each bar with the NaPS Atom - Sorbent Atom
            for j, value in enumerate(data[solvent]):
                ax.text(j + bar_width * i, value + 0.05, connections[solvent][j], ha='center', va='bottom', size=5)

        ax.set_xticks([p + bar_width for p in index])
        ax.set_xticklabels(x_tick_labels)
        y_min = min(min(values) for values in data.values())
        y_max = max(max(values) for values in data.values())
        y_range = y_max - y_min
        ax.set_ylim(1.5, y_max + 0.1 * y_range)
        ax.set_ylabel('Adsorption Distance (Å)', fontsize=15)
        ax.legend()
        ax.set_xlabel('NaPS', fontsize=15)
        
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
    x_tick_labels = [r"$Na_2S$", r"$Na_2S_2$", r"$Na_2S_4$", r"$Na_2S_6$", r"$Na_2S_8$"]
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
        fig, (ax_Na, ax_S) = plt.subplots(2, 1, figsize=(14, 10))
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
        ax_Na.set_xticklabels(x_tick_labels)
        ax_S.set_xticks([p + bar_width for p in index])
        ax_S.set_xticklabels(x_tick_labels)

        # set title for Na and S
        ax_Na.set_title('Na Oxidation States', fontsize=17)
        ax_S.set_title('S Oxidation States', fontsize=17)
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
        ax_Na.set_ylabel('Oxidation State', fontsize=15)
        ax_S.set_ylabel('Oxidation State', fontsize=15)
        ax_Na.legend(loc='best')
        ax_S.legend(loc='lower right')
        ax_Na.set_ylim(0.6, 1.2)
        if 'TiO2' in csv:
            ax_Na.set_ylim(0, 1.2)
        ax_Na.set_xlabel('NaPS', fontsize=15)
        ax_S.set_xlabel('NaPS', fontsize=15)
        # put gap between ax_Na and ax_S
        plt.subplots_adjust(hspace=0.3)
        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures/oxidation_states')
        os.makedirs(figures_dir, exist_ok=True)
        plt.savefig(os.path.join(figures_dir, f'{csv.split(".")[0]}.jpg'), bbox_inches='tight', pad_inches=0.1, dpi=dpi)

# --- This method plots the summary of binding energies for all materials --- #
def plot_binding_energy_summary():
    sorbents = ["graphene_vdw", "N_graphene", "FeN4_vdw", "TiO2", "NiS2", "NiS2_S"]
    # sorbents = ["graphene_vdw", "N_graphene", "FeN4_vdw", "NiS2", "NiS2_S"]
    sorbent_names = ["$Graphene$", "$N-graphene$", "$FeN_4-graphene$", "$TiO_2$", "$NiS_2$ Ni side", "$NiS_2$ S side"]
    # sorbent_names = ["$Graphene$", "$N-graphene$", "$FeN_4-graphene$", "$NiS_2$ Ni side", "$NiS_2$ S side"]
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/binding_energies')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv') and 'reference' not in csv]
    NaPS = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    NaPS_markers = ['o', 's', 'D', '^', 'v']
    solvents = ["Vacuum", "Glyme", "PC"]
    solvent_colors = ['b', 'r', 'g']
    # initialize dictionary
    data = {sorbent: {solvent: [] for solvent in solvents} for sorbent in sorbents}
    NaPS_marker_dict = {NaPS[i]: NaPS_markers[i] for i in range(len(NaPS))}
    solvent_color_dict = {solvents[i]: solvent_colors[i] for i in range(len(solvents))}
    # initialize plot
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    axs.set_xlabel('Sorbent')
    axs.set_ylabel('Binding Energy (kJ/mol)', fontsize=15)
    axs.set_xticklabels(sorbent_names)
    axs.set_xticks(range(len(sorbents)))
    for csv in csvs:
        sorbent = extract_sorbent(csv)
        if sorbent not in sorbents:
            continue
        csv_path = os.path.join(directory, csv)
        df = pd.read_csv(csv_path)
        # get sorbent index
        
        sorbent_index = sorbents.index(sorbent)
        for solvent in solvents:
            for NaPS_mol in NaPS:
                filtered_df = df[(df['NaPS'] == NaPS_mol) & (df['Solvent'] == solvent)]
                if not filtered_df.empty:
                    axs.plot(np.full(filtered_df.shape[0],sorbent_index), filtered_df['Binding Energy'].values[0]*conversion_factor, marker=NaPS_marker_dict[NaPS_mol], color=solvent_color_dict[solvent], markersize=10)
                else:
                    print(f"Missing data for {NaPS_mol} in {solvent} in {sorbent} in {csv}")

    
    # generate legend with section for NaPS and section for solvent
    NaPS_legend = [matplotlib.lines.Line2D([0], [0], marker=NaPS_markers[i], color='w', label=NaPS[i], markerfacecolor='black', markersize=10) for i in range(len(NaPS))]
    solvent_legend = [matplotlib.lines.Line2D([0], [0], marker='o', color='w', label=solvents[i], markerfacecolor=solvent_colors[i], markersize=10) for i in range(len(solvents))]

    # Create headers as text elements
    header_NaPS = matplotlib.lines.Line2D([0], [0], color='none', label='NaPS (marker)')
    header_solvent = matplotlib.lines.Line2D([0], [0], color='none', label='Solvent (color)')
    header_empty = matplotlib.lines.Line2D([0], [0], color='none', label='')

    # Combine headers and legend entries
    legend_elements = [header_NaPS] + NaPS_legend + [header_empty] + [header_solvent] + solvent_legend

    # Create the legend with a larger title
    legend = axs.legend(handles=legend_elements, loc='upper left', fontsize='x-large')

    plt.savefig('Figures/binding_energy_summary.jpg', dpi=dpi)


if __name__ == "__main__":
    # plot_cell_size_convergence()
    # plot_kpoint_convergence()
    plot_NaPS_energy()
    # plot_NaPS_oxidation()
    # plot_binding_energies()
    # plot_NaPS_bonds()
    # plot_adsorption_lengths()
    # plot_oxidation_states()
    # plot_binding_energy_summary()
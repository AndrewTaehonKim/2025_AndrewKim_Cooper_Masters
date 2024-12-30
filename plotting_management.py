import pandas as pd
import os

import matplotlib.pyplot as plt
from data_management import extract_sorbent

import matplotlib
matplotlib.use('Agg')  # Use a non-GUI backend

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
                    print(f"Missing data for {tick} in {solvent} in {csv}")
                    data[solvent].append(0)  # Handle missing data
        # Plot the data
        fig, ax = plt.subplots()
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
        figures_dir = os.path.join(os.getcwd(), 'Figures')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        plt.savefig(os.path.join(figures_dir, f'{extract_sorbent(csv)}.jpg'))
        plt.close()


# --- This method plots the bond lengths of the NaPS molecules in the adsorbed state vs. non-adsorbed state --- #
def plot_NaPS_bonds():
    # Define the directory containing the CSV files
    directory = os.path.join(os.getcwd(), 'Data-extracted/final/bonds')
    csvs = [csv for csv in os.listdir(directory) if csv.endswith('.csv')]
    # Define the x-ticks and solvents
    x_ticks = ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8"]
    solvents = ["Vacuum", "Glyme", "PC"]
    # Initialize a dictionary to hold the data
    # Read the data from the CSV files
    for csv in 

# --- This method plots the adsorption bond lengths between the NaPS molecules and the sorbent --- #
def plot_adsorption_lengths():
    pass


# --- This method plots the oxidation states of all molecules in the system and emphasizes the oxidation states of the NaPS and sites near the NaPS--- #
def plot_oxidation_states():
    pass


if __name__ == "__main__":
    # plot_binding_energies()
    plot_NaPS_bonds()
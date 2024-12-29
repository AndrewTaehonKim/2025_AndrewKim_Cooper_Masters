import pandas as pd
import os

import matplotlib.pyplot as plt
from data_management import extract_sorbent

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
                    data[solvent].append(-100)  # Handle missing data
        # Plot the data
        fig, ax = plt.subplots()
        bar_width = 0.2
        index = range(len(x_ticks))
        
        for i, solvent in enumerate(solvents):
            ax.bar([p + bar_width * i for p in index], data[solvent], bar_width, label=solvent, hatch=['/', '\\', ''][i])
        
        ax.set_xticks([p + bar_width for p in index])
        ax.set_xticklabels(x_ticks)
        ax.set_ylim(-0.005, 0.4)  
        ax.set_ylabel('Binding Energy (hartrees)')
        ax.legend()
        
        # Save the plot to a file
        figures_dir = os.path.join(os.getcwd(), 'Figures')
        if not os.path.exists(figures_dir):
            os.makedirs(figures_dir)
        plt.savefig(os.path.join(figures_dir, f'{extract_sorbent(csv)}.jpg'))
        plt.close()

if __name__ == "__main__":
    plot_binding_energies()
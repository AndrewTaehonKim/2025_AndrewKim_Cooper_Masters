import os

from data_management import *
from plotting_management import *

def extract_all_raw_data(run_all:bool=False, verbose:bool=False):
    # get working directory
    working_directory = os.getcwd()
    # Get list of categories
    categories = os.listdir(working_directory + "/Data-raw")
    if run_all:
        for category in categories:
            extract_data(category)
            if "@" in category:
                extract_adsorption_energies(category)
        print("All raw data extracted")
    else:
        csvs = [file for file in os.listdir(working_directory + "/Data-extracted/final/energies") if file.endswith('.csv')]
        for category in categories:
            if verbose:
                print(f"checking if {category} in {csvs}")
            if category+'.csv' not in csvs:
                extract_data(category)
                if verbose:
                    print(f"data extracted for {category}")
            else:
                if verbose:
                    print("skipped because found")
            if "@" in category:
                extract_adsorption_energies(category)

def plot_all_data():
    plot_binding_energies()
    print("All data plotted")
    return 0

# extract_all_raw_data()
extract_all_raw_data(run_all=True)

plot_all_data()

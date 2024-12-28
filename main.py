from data_management import *
import os

def extract_all_raw_data():
    # get working directory
    working_directory = os.getcwd()
    # Get list of categories
    categories = os.listdir(working_directory + "/Data-raw")
    for category in categories:
        extract_data(category)
    print("All raw data extracted")
extract_all_raw_data()
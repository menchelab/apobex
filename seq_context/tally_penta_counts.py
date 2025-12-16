import csv
import re
import argparse
import os

# Set up argument parsing
parser = argparse.ArgumentParser(description="Tally up sequence statistics.")
parser.add_argument("input_file", help="Path to the input file containing sequences.")
args = parser.parse_args()

# Input
input_file = args.input_file

# Extract sample ID from the file name
sample_ID = os.path.basename(input_file).split("_")[0]

output_file = f'{sample_ID}_sequence_context_statistics.csv'
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

# Initialize counters
conC = 0
conRTCA = 0
mutC = 0
mutRTCA = 0

# Define the patterns for conRTCA
#patterns = ["ATCA", "ATGA", "GTCA", "GTGA"]

# Define the patterns for conYTCA
patterns = ["CTCA", "CTGA", "TTCA", "TTGA"]

# Open and process the input file
with open(input_file, "r") as file:
    for line in file:
        # Skip header lines (lines starting with '>')
        if line.startswith(">"):
            continue
        
        # Remove newline characters and convert to uppercase
        sequence = line.strip().upper()
        
        # Count C and G bases
        conC += sequence.count("C") + sequence.count("G")
        
        # Count occurrences of the patterns
        for pattern in patterns:
            conRTCA += len(re.findall(pattern, sequence))

        # Check the 21st letter for mutC and mutRTCA
        if len(sequence) == 41:  # Ensure the sequence is 41 letters long
            context = sequence[18:22]  # 2 before + 1 after
            if context in ["ATCA", "ATGA", "GTCA", "GTGA"]:
                mutRTCA += 1
            else:
                mutC += 1

# Write results to a CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    
    # Write header
    writer.writerow(["Sample_ID","conC", "conRTCA", "mutC", "mutRTCA"])
    
    # Write results
    writer.writerow([sample_ID, conC, conRTCA, mutC, mutRTCA])

print(f"Results written to {output_file}")
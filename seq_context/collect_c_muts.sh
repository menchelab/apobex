#!/bin/bash

#SBATCH --output /nobackup/lab_menche/mmeyenberg/apobex/log_%j.log
#SBATCH --job-name=apobex
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.vcf"
    exit 1
fi

input_vcf="$1"

# Ensure the input VCF file exists
if [ ! -f "$input_vcf" ]; then
    echo "Error: Input file '$input_vcf' not found."
    exit 1
fi

# Generate the output file name
base_name=$(basename "$input_vcf" .vcf)
output_file="${base_name}_C_muts.txt"

# Process the VCF file
awk -v OFS="\t" '
    BEGIN { found_header = 0 }
    /^#CHROM/ { found_header = 1; next }
    found_header && $4 == "C" && ($5 == "G" || $5 == "T") {
        print $1, $2, $4, $5
    }
' "$input_vcf" > "$output_file"

echo "Filtered data written to '$output_file'."
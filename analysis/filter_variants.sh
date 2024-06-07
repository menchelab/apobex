#!/bin/bash

##########################################################################################
# This script takes in a folder path to the directory containing the vcf files           #
# to be filtered and a path to the vcf containing the variants found in                  #
# the WGS of the parental clone.                                                         #   
#                                                                                        #
# It sequentially applies the following filters and automatically deletes                #
# intermedidiary files:                                                                  #   
#                                                                                        #   
# Filter 1 = Keep only variants with Mutect2 filter status == PASS                       #
# Filter 2 = Keep only variants with a VAF >= 0.3 and <= 0.7                             #
# Filter 3 = Keep only variants that were mapped to canonical chromosomes (1-22, X/Y)    #
# Filter 4 = Discard varients which are also present in the WGS of the parental clone    #
#                                                                                        #
# Output = new directory containing vcf files with all filters applied                   #
##########################################################################################


# Prompt user for input
read -p "Enter path to raw vcfs: " sample_path
read -p "Enter path to WGS vcf of parental clone: " parental_path

# Create output directory

mkdir "${sample_path}"/filtered_vcf
out_path="${sample_path}"/filtered_vcf

# Apply Filter 1

cd $sample_path

#---------------------------------------------
# Apply Filter 1

echo "Applying Filter 1"

for file in "${sample_path}"/*.vcf
    do
        # Extract the base name without extension
        base_name=$(basename "$file" .vcf)
        
        # Use awk to filter lines where filter status == PASS and save to a new file
        awk '$1 ~ /^#/ || $7 == "PASS" ' "$file" > "${out_path}/${base_name}_pass.vcf"
    done


echo "Variants with filter status PASS retained"

#---------------------------------------------
# Apply Filter 2

echo "Applying Filer 2"

for file in "${out_path}"/*.vcf
    do
        # Extract the base name without extension
        base_name=$(basename "$file" .vcf)

        # Filter variants based on VAF using bcftools view
         bcftools view -i "AF >= 0.3 & AF <= 0.7" -o "$out_path/${base_name}_VAF.vcf" "$file"
    done

# Clean output directory
cd ${out_path}
rm *_pass.vcf

echo "Retained variants with VAF >= 0.3 and <= 0.7"

#---------------------------------------------
# Apply Filter 3

echo "Applying Filter 3"

# Prepare file input
for file in "${out_path}"/*.vcf
    do 
        bgzip $file
    done

for file in "${out_path}"/*.vcf.gz
    do
        bcftools index $file
    done

# Filter by Chromsomes
chromosomes=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,\
chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY

for file in "${out_path}"/*.vcf.gz
    do
        # Define output file name
        base_name=$(basename "$file" .vcf)
        output_vcf="$out_path/${base_name}_chrom.vcf"

        # Filter VCF by Chromosomes
        bcftools view $file --regions $chromosomes -o $output_vcf
    done

# Clean out_path
cd ${out_path}
rm *.vcf.gz
rm *.vcf.gz.csi

echo "Retained variants on canonical chromsomes"

#---------------------------------------------
# Apply Filter 4

echo "Applying Filter 4"

for file in "${out_path}"/*.vcf
    do  
        # Extract the base name without extension
        base_name=$(basename "$file" .vcf)

        # Prep output file
        sed -n '/^#CHROM/q;p' $file > $out_path/${base_name}_WGS.vcf
        grep "^#CHROM" $file >> $out_path/${base_name}_WGS.vcf

        # Apply WGS subtractionls
        bedtools subtract -A -a $file -b ${parental_path} >> $out_path/${base_name}_WGS.vcf

    done

# Clean out_path
cd ${out_path}
rm *_chrom.vcf

# Simplify file names
for file in "${out_path}"/S*_*.vcf
    do
        base_name=$(basename "$file")
        sample_id=$(echo "$base_name" | awk -F'_' '{print $1}')
        
        # Construct the new file name
        new_name="${sample_id}_filtered.vcf"
        
        # Rename the file
        mv "$file" "$new_name"
    done


echo "Eliminated Variants which are also present in the parental clone"
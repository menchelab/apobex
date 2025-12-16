library(readxl)
library(vcfR)
library(stringr)
library(dplyr)
install.packages("tidyr")
library(tidyr)

#-------------------------------------------------------------
####INPUT####

#read PCAWG sample sheet and filter for "WGS" and "Whitelist"
id_map <- read.table("./pcawg_sample_sheet.tsv", sep = "\t", header = TRUE)
id_map <- id_map[id_map$library_strategy=="WGS",]
id_map <- id_map[id_map$donor_wgs_exclusion_white_gray =="Whitelist",]

#read "PCAWG_sigProfiler_SBS_signatures_in_samples"
samples <- read.csv("./PCAWG_sigProfiler_SBS_signatures_in_samples.csv", sep=";")

#merge the two tables to combine sample and mutational signature information
samples_m <- merge(id_map, samples, by.x=c("icgc_specimen_id"), by.y=c("Sample.Names"))
samples_m <- samples_m[samples_m$library_strategy=="WGS",]
samples_m <- samples_m[samples_m$donor_wgs_exclusion_white_gray =="Whitelist",]


#-------------------------------------------------------------
####COUNT NUMBER OF MUTATIONS PER SAMPLE####

#define function to count total number of mutations per sample
count_mutations <- function(samples_mapped){
  
  samples_mapped[, "mutation_counts"] <- NA
  
  #identify snvs in genes
  vcf_list_snv <- list.files("./PCAWG_combined_data/all_vcf_snv_extracted/")
  counter <- 0
  total <- 0
  
  for(i in 1:nrow(samples_mapped)) {
    row <- samples_mapped[i,]
    id <- row$aliquot_id
    
    if(total%%100 == 0){
      print(total)
    }
    
    for(j in 1:length(vcf_list_snv)) {
      filename <- vcf_list_snv[j]
      if(startsWith(filename, id)){
        vcf <- read.vcfR(paste("./PCAWG_combined_data/all_vcf_snv_extracted/", filename, sep = ""), verbose = FALSE)
        vcf_df <- as.data.frame(vcf@fix)
        vcf_df <- transform(vcf_df, POS = as.numeric(POS))
        vcf_df <- vcf_df[vcf_df$FILTER == ".", ]
        
        if(nrow(vcf_df) == 0){
          counter <- counter + 1
          row["mutation_counts"] = FALSE
        } else {
          row["mutation_counts"] = nrow(vcf_df)
        }
        
        total <- total + 1
        break
      }
    }
    samples_mapped[i,] <- row
  }
  
  print(vcf_df$INFO)
  return(samples_mapped)
}

#execute function
samples_m_counted <- count_mutations(samples_m)

#-------------------------------------------------------------
####IDENTIFY SNVs & INDELs IN GENES OF INTEREST####

#define function to identify snv's and indels for a list of genes 
find_mutant <- function(samples_mapped, gene_name, gene_dict){
 
  samples_mapped[,paste(gene_name, "_mut")] <- NA
  samples_mapped[,paste(gene_name, "_indel")] <- NA
  
  r <- gene_dict[gene_dict$Name==gene_name,]
  Chrom = r[,"Chrom"]
  Start = r[,"Start"]
  End = r[,"End"]
  

  #identify snvs in genes
  vcf_list_snv <- list.files("./all_vcf_snv_extracted/")
  counter <- 0
  total <- 0
  for(i in 1:nrow(samples_mapped)) {
    row <- samples_mapped[i,]
    id <- row$aliquot_id
    if(total%%100==0){
      print(total)
    }
    
    for(j in 1:length(vcf_list_snv)) {
      filename <- vcf_list_snv[j]
      if(startsWith(filename, id)){
        vcf <- read.vcfR(paste("./all_vcf_snv_extracted/",filename, sep= ""), verbose = FALSE )
        vcf_df <- as.data.frame(vcf@fix)
        vcf_df <- transform(vcf_df, POS = as.numeric(POS))
        vcf_df$Variant_Class <- str_split_fixed(vcf_df$INFO, ";Variant_Classification=", 2)
        vcf_df <- vcf_df[vcf_df$CHROM==Chrom,]
        vcf_df <- vcf_df[vcf_df$POS>=Start,]
        vcf_df <- vcf_df[vcf_df$POS<=End,]
        vcf_df <- vcf_df[vcf_df$Variant_Class[,2]!="Intron",]
        vcf_df <- vcf_df[vcf_df$FILTER==".",]
        
        if(nrow(vcf_df) == 0){
          counter <- counter + 1
          row[,paste(gene_name, "_mut")] = FALSE
        }
        else{
          row[,paste(gene_name, "_mut")] = TRUE
        }
        
        
        total <- total + 1
        break
      }
    }
    samples_mapped[i,] = row
    
  }
  
  #identify indels in genes
  vcf_list <- list.files("./all_vcf_indel_extracted/")
  counter <- 0
  total <- 0
  for(i in 1:nrow(samples_mapped)) {
    row <- samples_mapped[i,]
    id <- row$aliquot_id
    if(total%%100==0){
      print(total)
    }
    
    for(j in 1:length(vcf_list)) {
      filename <- vcf_list[j]
      if(startsWith(filename, id)){
        vcf <- read.vcfR(paste("./all_vcf_indel_extracted/",filename, sep= ""), verbose = FALSE )
        vcf_df <- as.data.frame(vcf@fix)
        vcf_df <- transform(vcf_df, POS = as.numeric(POS))
        vcf_df$Variant_Class <- str_split_fixed(vcf_df$INFO, ";Variant_Classification=", 2)
        vcf_df <- vcf_df[vcf_df$CHROM==Chrom,]
        vcf_df <- vcf_df[vcf_df$POS>=Start,]
        vcf_df <- vcf_df[vcf_df$POS<=End,]
        vcf_df <- vcf_df[vcf_df$Variant_Class[,2]!="Intron",]
        vcf_df <- vcf_df[vcf_df$FILTER==".",]
        
        if(nrow(vcf_df) == 0){
          counter <- counter + 1
          row[,paste(gene_name, "_indel")] = FALSE
        }
        else{
          row[,paste(gene_name, "_indel")] = TRUE
        }
        
        
        total <- total + 1
        break
      }
    }
    samples_mapped[i,] = row
    
  }
  print(vcf_df$INFO)
  return(samples_mapped)
  
}

#execute the funcion for a list of genes defined in gene_dict
#define genes of interest and their genomic location in gene_dict; load the file
gene_dict <- read.csv("./gene_dict.csv", sep=";")

for(k in 1:nrow(gene_dict)){
  gene_name <- gene_dict[k,"Name"]
  print(gene_name)
  samples_m_counted <- find_mutant(samples_m_counted, gene_name, gene_dict)
}

#remove whitespace in column names
colnames(samples_m_counted) <- gsub("[[:space:]]", "", colnames(samples_m_counted))


#-------------------------------------------------------------
####ANALYSIS: define wt and alt (either mut or indel = TRUE)####

#rename data table
samples_genes <- samples_m_counted

#remove whitespace from column names
colnames(samples_genes) <- gsub("[[:space:]]", "", colnames(samples_genes))

#specify if gene is wt or not in new columns, remove columns with NA
for (l in 1:nrow(gene_dict)) {
  gene_name <- gsub("[[:space:]]", "", gene_dict[l, "Name"])  # Remove all whitespace
  
  # Construct column names (without space)
  mut_col <- paste(gene_name, "_mut", sep = "")
  indel_col <- paste(gene_name, "_indel", sep = "")
  wt_col <- paste(gene_name, "_wt", sep = "")
  alt_col <- paste(gene_name, "_alt", sep = "")
  
  # Further specification
  samples_genes[[wt_col]] <- samples_genes[[mut_col]] == FALSE & samples_genes[[indel_col]] == FALSE
  samples_genes[[alt_col]] <- samples_genes[[mut_col]] == TRUE | samples_genes[[indel_col]] == TRUE
  
  # Remove rows with NA
  samples_genes[[wt_col]] <- samples_genes[[wt_col]] & !is.na(samples_genes[[wt_col]])
  samples_genes[[alt_col]] <- samples_genes[[alt_col]] & !is.na(samples_genes[[alt_col]])
}

#check if column names have whitespace
print(colnames(samples_genes))


#create the new columns "E3_alt" and "E3_wt" based on the specified conditions
samples_genes$E3_alt <- samples_genes$UBR4_alt | samples_genes$UBR5_alt | samples_genes$HUWE1_alt
samples_genes$E3_wt <- samples_genes$UBR4_wt & samples_genes$UBR5_wt & samples_genes$HUWE1_wt

#filter the data for E3_wt samples
E3_wt_data <- samples_genes[samples_genes$E3_wt == TRUE, ]

#add genotype column for E3_wt samples
E3_wt_data$genotype <- "E3_wt"

#filter the data for E3_alt samples
E3_alt_data <- samples_genes[samples_genes$E3_alt == TRUE, ]

#add genotype column for E3_alt samples
E3_alt_data$genotype <- "E3_alt"

#assign the data frames to the R environment
assign("samples_genes_E3_wt", E3_wt_data)
assign("samples_genes_E3_alt", E3_alt_data)


#-------------------------------------------------------------
####ANALYSIS####
##prepare data to visualize##
#iterate over each gene name from gene_dict
for (i in 1:nrow(gene_dict)) {
  # Access the gene name
  gene_name <- gene_dict[i, "Name"]
  
  # Define the column names for wt and alt samples of the current gene
  wt_col <- paste0(gene_name, "_wt")
  alt_col <- paste0(gene_name, "_alt")
  
  # Filter the data for wt and alt samples of the current gene
  wt_data <- samples_genes[samples_genes[[wt_col]] == TRUE, ]
  alt_data <- samples_genes[samples_genes[[alt_col]] == TRUE, ]
  
  # Add genotype column to wt and alt data frames
  wt_data$genotype <- paste(gene_name, "wt", sep = "_")
  alt_data$genotype <- paste(gene_name, "alt", sep = "_")
  
  # Assign the data frames to the R environment
  assign(paste0("samples_genes_", gene_name, "_wt"), wt_data)
  assign(paste0("samples_genes_", gene_name, "_alt"), alt_data)
}

#filter the data for E3_wt samples
E3_wt_data <- samples_genes[samples_genes$E3_wt == TRUE, ]
#add genotype column for E3_wt samples
E3_wt_data$genotype <- "E3_wt"

#filter the data for E3_alt samples
E3_alt_data <- samples_genes[samples_genes$E3_alt == TRUE, ]
#add genotype column for E3_alt samples
E3_alt_data$genotype <- "E3_alt"

#assign the data frames to the R environment
assign("samples_genes_E3_wt", E3_wt_data)
assign("samples_genes_E3_alt", E3_alt_data)



#-------------------------------------------------------------
####ANALYSIS: normalization to total number of mutations in sample####

#remove whitespace from column names
colnames(samples_genes) <- gsub("[[:space:]]", "", colnames(samples_genes))

#list of columns starting with "SBS"
sbs_columns <- colnames(samples_genes)[grepl("^SBS", colnames(samples_genes))]

#create a copy of samples_genes for normalization
samples_genes_normalized <- samples_genes

#normalize each SBS column to the total number of mutations
for (col in sbs_columns) {
  samples_genes_normalized[[col]] <- samples_genes_normalized[[col]] / samples_genes_normalized$mutation_counts
}


#-------------------------------------------------------------
####ANALYSIS: prepare data to visualize - WITH NORMALIZATION#### 

#iterate over each gene name from gene_dict
for (i in 1:nrow(gene_dict)) {
  # Access the gene name
  gene_name <- gene_dict[i, "Name"]
  
  # Define the column names for wt and alt samples of the current gene
  wt_col <- paste0(gene_name, "_wt")
  alt_col <- paste0(gene_name, "_alt")
  
  # Filter the data for wt and alt samples of the current gene
  wt_data <- samples_genes_normalized[samples_genes_normalized[[wt_col]] == TRUE, ]
  alt_data <- samples_genes_normalized[samples_genes_normalized[[alt_col]] == TRUE, ]
  
  # Add genotype column to wt and alt data frames
  wt_data$genotype <- paste(gene_name, "wt", sep = "_")
  alt_data$genotype <- paste(gene_name, "alt", sep = "_")
  
  # Assign the data frames to the R environment
  assign(paste0("samples_genes_normalized_", gene_name, "_wt"), wt_data)
  assign(paste0("samples_genes_normalized_", gene_name, "_alt"), alt_data)
}

#filter the data for E3_wt samples
E3_wt_data <- samples_genes_normalized[samples_genes_normalized$E3_wt == TRUE, ]
#add genotype column for E3_wt samples
E3_wt_data$genotype <- "E3_wt"

#filter the data for E3_alt samples
E3_alt_data <- samples_genes_normalized[samples_genes_normalized$E3_alt == TRUE, ]
#add genotype column for E3_alt samples
E3_alt_data$genotype <- "E3_alt"

#assign the data frames to the R environment
assign("samples_genes_normalized_E3_wt", E3_wt_data)
assign("samples_genes_normalized_E3_alt", E3_alt_data)

#manually define which genes to plot
analysis_df <- rbind(samples_genes_normalized_UBR4_wt,
                     samples_genes_normalized_UBR4_alt,
                     samples_genes_normalized_UBR5_wt,
                     samples_genes_normalized_UBR5_alt,
                     samples_genes_normalized_HUWE1_wt,
                     samples_genes_normalized_HUWE1_alt,
                     samples_genes_normalized_E3_wt,
                     samples_genes_normalized_E3_alt)



analysis_df <- rbind(samples_genes_normalized_HECTD4_wt,
                     samples_genes_normalized_HECTD4_alt,
                     samples_genes_normalized_NEDD4L_wt,
                     samples_genes_normalized_NEDD4L_alt,
                     samples_genes_normalized_HECTD1_wt,
                     samples_genes_normalized_HECTD1_alt,
                     samples_genes_normalized_PCGF3_wt,
                     samples_genes_normalized_PCGF3_alt)




analysis_df <- rbind(samples_genes_normalized_UBR4_wt,
                     samples_genes_normalized_UBR4_alt,
                     samples_genes_normalized_UBR5_wt,
                     samples_genes_normalized_UBR5_alt,
                     samples_genes_normalized_HUWE1_wt,
                     samples_genes_normalized_HUWE1_alt,
                     samples_genes_normalized_E3_wt,
                     samples_genes_normalized_E3_alt,
                     samples_genes_normalized_HECTD4_wt,
                     samples_genes_normalized_HECTD4_alt,
                     samples_genes_normalized_NEDD4L_wt,
                     samples_genes_normalized_NEDD4L_alt,
                     samples_genes_normalized_HECTD1_wt,
                     samples_genes_normalized_HECTD1_alt,
                     samples_genes_normalized_PCGF3_wt,
                     samples_genes_normalized_PCGF3_alt)


#-------------------------------------------------------------
####ANALYSIS: Visualize data####
library(ggplot2)
library(ggpubr)


#plot mutational_counts or signature of interest
SBS <- "SBS2"

#manually define which genes to plot
p <- ggboxplot(analysis_df, x = "genotype", y = SBS, color = "genotype",  add = "jitter")+
  scale_y_continuous(trans = 'log10') +
  scale_color_manual(values = c("UBR4_wt" = "#71A6DA", "UBR4_alt" = "#16365B", 
                                "UBR5_wt" = "#71A6DA", "UBR5_alt" = "#16365B",
                                "HUWE1_wt" = "#71A6DA", "HUWE1_alt" = "#16365B",
                                "E3_wt" = "#71A6DA", "E3_alt" = "#16365B",
                                "HECTD4_wt" = "#71A6DA", "HECTD4_alt" = "#16365B", 
                                "NEDD4L_wt" = "#71A6DA", "NEDD4L_alt" = "#16365B",
                                "HECTD1_wt" = "#71A6DA", "HECTD1_alt" = "#16365B",
                                "PCGF3_wt" = "#71A6DA", "PCGF3_alt" = "#16365B"))   # Manual color specification for each genotype

#show the plot
p


#assuming you have already performed the comparison of means
wilcox <- compare_means(formula = SBS2 ~ genotype, data = analysis_df, method = "wilcox.test", paired = FALSE,
                        group.by = NULL, ref.group = NULL)


print(wilcox, n=120)

#filter significant comparisons (I don't care about ns)
significant_comparisons <- wilcox[wilcox$p.signif != "ns", ]

#print significant comparisons
print(significant_comparisons, n = 1000)

##only show comparisons between wt and alt for the same gene##
#extract gene names from group1 and group2
gene_name1 <- str_extract(significant_comparisons$group1, "^[^_]+")
gene_name2 <- str_extract(significant_comparisons$group2, "^[^_]+")

#filter significant comparisons between wt and alt from the same name
significant_comparisons_same_name <- significant_comparisons[
  gene_name1 == gene_name2, 
]

#print significant comparisons
print(significant_comparisons_same_name, n = 1000)

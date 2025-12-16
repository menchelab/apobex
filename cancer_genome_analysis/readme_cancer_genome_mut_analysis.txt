Analysis of cancer genome sequencing data

This analysis tests the correlation between E3 ligase mutations and APOBEC signature burden in cancer samples.
Whole-genome sequencing (WGS) samples used in the PCAWG study, as well as the SBS mutational signature analysis from the same study were used.
Cancer samples, in which the E3 ligase of interest is mutated (at least one SNV or InDel outside intronic regions) were identified.
SBS mutational signatures were normalized to the total number of signatures in the respective sample.
The correlation between the E3 ligase genotype and the SBS mutational burden was plotted and analyzed.


Public and restricted cancer whole-genome sequencing data from ICGC and TCGA datasets were obtained.
Public data (vcf) were obtained through the ICGC Data Portal (https://dcc.icgc.org/releases/PCAWG): DCC/PCAWG/consensus_snv_indel/final_consensus_snv_indel_passonly_icgc.public.tgz
Restricted data (vcf) was obtained through the NCBI dbGaP with granted access by the NIH.
SBS mutational signatures analyzed in the PCAWG study were obtained from the ICGC Data Portal: DCC/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv
The PCAWG sample sheet was obtained from the ICGC Data Portal: /PCAWG/data_releases/latest/pcawg_sample_sheet.tsv



Session info:

R version 4.1.3 (2022-03-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] LC_COLLATE=English_Austria.1252  LC_CTYPE=English_Austria.1252    LC_MONETARY=English_Austria.1252 LC_NUMERIC=C                    
[5] LC_TIME=English_Austria.1252    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tidyr_1.3.0          dplyr_1.1.2          stringr_1.5.1        vcfR_1.14.0          readxl_1.4.2         ggpubr_0.6.0         ggplot2_3.4.4       
 [8] GenomicRanges_1.46.1 GenomeInfoDb_1.30.1  IRanges_2.28.0       S4Vectors_0.32.4     BiocGenerics_0.40.0 

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0       purrr_1.0.1            splines_4.1.3          lattice_0.21-8         carData_3.0-5          colorspace_2.1-0      
 [7] vctrs_0.6.1            generics_0.1.3         viridisLite_0.4.2      mgcv_1.8-42            utf8_1.2.3             rlang_1.1.0           
[13] pillar_1.9.0           glue_1.6.2             withr_3.0.0            GenomeInfoDbData_1.2.7 lifecycle_1.0.4        zlibbioc_1.40.0       
[19] munsell_0.5.0          ggsignif_0.6.4         gtable_0.3.4           cellranger_1.1.0       permute_0.9-7          parallel_4.1.3        
[25] fansi_1.0.4            broom_1.0.5            Rcpp_1.0.10            pinfsc50_1.3.0         backports_1.4.1        scales_1.2.1          
[31] vegan_2.6-4            XVector_0.34.0         abind_1.4-5            farver_2.1.1           digest_0.6.31          stringi_1.7.12        
[37] rstatix_0.7.2          grid_4.1.3             cli_3.6.1              tools_4.1.3            bitops_1.0-7           magrittr_2.0.3        
[43] RCurl_1.98-1.12        tibble_3.2.1           cluster_2.1.4          ape_5.7-1              car_3.1-2              pkgconfig_2.0.3       
[49] Matrix_1.5-4           MASS_7.3-58.3          rstudioapi_0.15.0      R6_2.5.1               nlme_3.1-162           compiler_4.1.3
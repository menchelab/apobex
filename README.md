# apobex

This repository contains all data and code to analyze and visualize the mutational data presented in:
**"Ubiquitin E3 ligases UBR4, UBR5 and HUWE1 protect human genome integrity by targeting cancer-associated APOBEC3 deaminases for degradation"**
link to preprint: https://www.biorxiv.org/content/biorxiv/early/2024/04/26/2024.04.23.590688.full.pdf

## Installation instructions

#### 1. clone the GitHub repository:
```bash
git clone https://github.com/menchelab/apobex.git
```
#### 2. Navigate into the ApobeX directory within the cloned repository
```bash
cd apobex #go into cloned repository
```
#### 3. Install the ApobeX package from the local repository:
```R
library(devtools)
install("ApobeX")
```
#### 4. Load the ApobeX package:
```R
library(ApobeX)
```

### Analysis
run notbook ../analysis/bootstrapped_background_correction.ipynb

### Plotting

#### figures
run notebook ../plotting/notbooks/3_apobex_plotting_refitting.ipynb

#### signature presence test
run norebook ../plotting/notbooks/3_apobex_msigact.ipynb

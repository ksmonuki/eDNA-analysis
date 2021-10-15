# eDNA Captures Depth Partitioning in a Kelp Forest Ecosystem

Pre-print is available [here](https://www.biorxiv.org/content/10.1101/2021.06.01.446542v1)

Keira Monuki*, Paul H Barber, Zachary Gold

Affiliation: University of California, Los Angeles, Los Angeles, CA

*Corresponding Author (Current address is Bodega Marine Laboratory, Bodega Bay, CA)

### Abstract 

Environmental DNA (eDNA) metabarcoding is an increasingly important tool for surveying biodiversity in marine ecosystems. However, the scale of temporal and spatial variability in eDNA signatures, and how this variation may impact eDNA-based marine biodiversity assessments, remains uncertain. To address this question, we systematically examined variation in vertebrate eDNA signatures across depth (0 m to 10 m) and horizontal space (nearshore kelp forest and surf zone) over three successive days in Southern California. Across a broad range of teleost fish and elasmobranchs, results showed significant variation in species richness and community assemblages between surface and depth, reflecting microhabitat depth preferences of common Southern California nearshore rocky reef taxa. Community assemblages between nearshore and surf zone sampling stations at the same depth also differed significantly, consistent with known habitat preferences. Additionally, assemblages also varied across three sampling days, but 69% of habitat preferences remained consistent. Results highlight the sensitivity of eDNA in capturing fine-scale vertical, horizontal, and temporal variation in marine vertebrate communities, demonstrating the ability of eDNA to capture a highly localized snapshot of marine biodiversity in dynamic coastal environments.

### Description

This page is dedicated to hosting the code for the accepted manuscript. 

Included on this page is: 

#### 1. Scripts used to conduct analyses
i. *decontam/keira_decontamination_updated_080520.Rmd* Code to decontaminate the taxonomic tables 
ii. *decontam/SOM_keira_08062020.R* Code to run the site occupancy model
iii. *eDNA_analysis_code.Rmd* Code to run downstream alphadiversity and betadiversity analyses

#### 2. Data

The data is available on [Dryad](https://doi.org/10.5068/D18H47)

### Dependencies 

Here are the main dependencies required. The additional R packages not listed below are included in the *eDNA_analysis_code.Rmd*

- R (code was run using R version 4.0.3)
- RStudio
- The following packages: 
  * tidyverse
  * rstan
  * shinystan
  * bayesplot
  * broom
  * vegan
  * proxy
  * reshape2
  * microDecon
  * gradientForest
  * extendedForest
  * ggplot2
  * phyloseq
  * ranacapa

### Workflow

1. The first script to run is *decontam/keira_decontamination_updated_080520.Rmd*. The input files needed for this script are *fishcard_12S_all_ASV_raw_taxonomy_60_edited.txt*, *KM_metadata_4.2.19_edited.txt*, *best_hashes_chosen_08052020*. 

    Once you complete the MicroDecon Miu step and create the output file *Cleaning.before.Occ.model*, move to the second step. 

2. The next script to run is the site occupancy model *decontam/SOM_keira_08062020.R*. The input files needed for this script are *Cleaning.before.Occ.model*, *KM_metadata_4.2.19_edited.txt*, *best_hashes_chosen_08052020*. 

    Once you run the script and obtain output files, move back to the decontamination script from Step 1. 

3. In the script *decontam/keira_decontamination_updated_080520.Rmd*, move to Cleaning Process 4 (Line 407). You will need two output files created from the site occupancy script in Step 2: *keira_occupancy_results_all.RDS* and *keira_occupancy_results_9.RDS*. Finish running the script. 

4. To run the analyses, use code *eDNA_analysis_code.Rmd*. You will need the input files *KM_metadata_8.7.20.txt*, *post_occupancy_results_merged_runs_separate_eDNA_index.RDS* (output from decontamination script in decontam folder), *ASV.summary_post_occ.RDS* (output from decontamination script in decontam folder), *KM_metadata_gradientForest_offshore_4.1.20.txt*, *KM_metadata_gradientForest_2.13.20.txt*, *ASV_post_occ_sum_taxonomy_pres_abs_separate.csv* (output from decontamination script in decontam folder). 


#### If you have any questions, please reach out to Keira Monuki at ksmonuki@ucdavis.edu










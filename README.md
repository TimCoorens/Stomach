## Code and data repo
This repository contains the code and filtered data accompanying the paper "The somatic mutation landscape of normal gastric epithelium", currently available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.03.17.585238v1).
  

For any queries, raise an issue on the GitHub page or email Tim Coorens (tcoorens@broadinstitute.org).

### Data
The folder Sequoia_output contains the filtered depth and variant depth output matrices from Sequoia ("*_NR_filtered_all.txt" and "*_NV_filtered_all.txt", respectively), tree files ("*_tree_with_branch_length.tree") and mutation mapping to trees ("*_assigned_to_branches.txt") for the WGS data of all donors, for both somatic SNVs and indels. Other data necessary for reproducibility can be found in the copy of the Extended Data Tables (Extended_Data_Table_1_7_R1_final.xlsx).

The folder HDP_signatures contains the extracted indel and SNV signatures.

### Scripts
"Mutation_burden_analyses.R" contains the code to generate Figure 1, Extended Data Figure 1 and Figure 2, along with all statistical tests in the section "Mutation rates of normal gastric epithelium". 

"Signature_analyses.R" contains the code to generate Figure 3 and Extended Data Figures 5-6, along with all statistical tests in the section "Mutational signatures and processes in normal gastric epithelium". 

"CNV_SV_analyses.R" contains the code to generate Figure 4 and Extended Data Figures 7-8, along with all statistical tests in the section "Recurrent trisomies in normal gastric glands".

"Driver_analyses.R" contains the code to generate Figure 5 and Extended Data Figure 9, along with all statistical tests in the section "Driver mutations in normal gastric glands"

### System requirements
All code has been run using R (v4.3.1) and will run on the provided input files. Note that R packages loaded in individual scripts need to installed prior to use. 

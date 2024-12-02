## Code and data repo

The R script in this repository takes in matrices of variant depth and total depth across from multiple clonal samples of the same donor, filters out germline and artefactual variants, constructs the phylogenetic tree topology and maps somatic mutations to the tree branches using the [TreeMut](https://github.com/nangalialab/treemut) package. 
For any queries, raise an issue on the GitHub page or email Tim Coorens (tcoorens@broadinstitute.org) or Mike Spencer Chapman (ms56@sanger.ac.uk).

### Data
The folder Sequoia_output contains the filtered NR and NV output matrices from Sequoia and tree files for the WGS data of all donors, for both somatic SNVs and indels. Other data necessary for reproducibility can be found in the copy of the Extended Data Tables (Extended_Data_Table_1_7_R1_final.xlsx).

### Scripts
"Mutation_burden_analyses.R" contains the code to generate Figure 1, Extended Data Figure 1 and Figure 2, along with all statistical tests in the section "Mutation rates of normal gastric epithelium". 
"Signature_analyses.R" contains the code to generate Figure 3 and Extended Data Figures 5-6, along with all statistical tests in the section "Mutational signatures and processes in normal gastric epithelium". 
"CNV_SV_analyses.R" contains the code to generate Figure 4 and Extended Data Figure 7, along with all statistical tests in the section "Recurrent trisomies in normal gastric glands".
"Driver_analyses.R" contains the code to generate Figure 5 and Extended Data Figures 8, along with all statistical tests in the section "Driver mutations in normal gastric glands"
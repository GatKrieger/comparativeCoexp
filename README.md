# Comparison of co-expression networks between species
code for co-expression analysis for the paper "Independent evolution of transcript abundance and regulatory dynamics". 

# Pipeline
The script "compare_PCMs.m" is a pipeline that runs the other functions, coded in matlab.
# Input file
The input is a table of normalized read counts in excel, where each sheet contains gene expression samples collected from a single species.
# Test input
"expression_per_species_test.xlsx"
This is a subset of the data presented in the article.
Full data will be avilable soon.
# Output files
The similairity in co-expression pattern of a gene between species is termed regulatory similarity.
"divT_<date>.mat" is a matlab table that lists the regulatory similarity scores of each gene in the desired comparisons.
"ggc_<date>.mat" is a matlab struct that contains the pairwise correlation matrices (PCMs) of each species.
# Pre-requisits
Matlab
# Citing
Independent evolution of transcript abundance and gene regulatory dynamics
Gat Krieger, Offir Lupo, Avraham A. Levy, Naama Barkai
bioRxiv 2020.01.22.915033; doi: https://doi.org/10.1101/2020.01.22.915033


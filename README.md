# GenomeMUSter
GenomeMUSter is available at https://muster.jax.org/ 

Please cite:
GenomeMUSter mouse genetic variation service enables multi-trait, multi-population data integration and analyses
Robyn L. Ball*, Molly A. Bogue, Hongping Liang, Anuj Srivastava, David G. Ashbrook, Anna Lamoureux, Matthew W. Gerring, Alexander S. Hatoum, Matthew Kim, Hao He, Jake Emerson, Alexander K. Berger, David O. Walton, Keith Sheppard, Baha El Kassaby, Francisco Castellanos, Govind Kunde-Ramamoorthy, Lu Lu, John Bluis, Sejal Desai, Beth A. Sundberg, Gary Peltz, Zhuoqing Fang, Gary A. Churchill, Robert W. Williams, Arpana Agrawal, Carol J. Bult, Vivek M. Philip, Elissa J. Chesler.
bioRxiv 2023.08.08.552506; doi: https://doi.org/10.1101/2023.08.08.552506

*Robyn L. Ball, PhD (email: robyn.ball at jax dot org)

Data from sixteen variation datasets (Table 1)  were merged by chromosome and within chromosome in 10 Mb segmented regions. Once the datasets were merged and harmonized, missing genotypes were imputed for each strain in each 10 Mb region.

Code is provided for a segment of the genome, i.e. chr, start, end. For purposes of demonstration, we provide data for a small region.

Use **MUSter_example.R** for a walkthrough of the process and refer to notes below.

## MERGE DATASETS

This pipeline merges across datasets for a specified genomic region (chr, start, end). Genotype data are expected to be in .csv format with chr, bp38 taken as the genomic coordinates followed by strain names as column names s.t. each row corresponds to a unique chr, bp38 site with genotypes A, C, T, G, H, or "".

This pipeline merges across datasets for a specified genomic region (chr, start, end). Genotype data are expected to be in .csv format with chr, bp38 taken as the genomic coordinates followed by strain names as column names s.t. each row corresponds to a unique chr, bp38 site with genotypes A, C, T, G, H, or "".

1. run merge_datasets.R

The output for chr, start, end is as follows:
- merged_chr{chr}_{start}-{end}.csv, e.g., merged_chr1_103000029-113000029.csv contains the merged dataset across all 11 datasets for the genomic region
- imputed_prob_chr{chr}_{start}-{end}.csv is the same shape as merged_chr{chr}_{start}-{end}.csv and contains either 1 (genotype is in the dataset) or 0 (genotype is not in the dataset). This file will be updated in the imputation step. If the genotype is imputed it will be between 0 and 1 (0,1) such that 1 indicated a typed/known SNP.
- mismatches/mismatches_chr{chr}_{start}-{end}_{dataset}.csv, e.g., mismatches_chr1_103000029-113000029_Broad2.csv, contains the genotype call for a strain and genomic location in the dataset if it was different than the genotype call in the merged dataset. These mismatched calls are collected as the datasets are merged so we can use them to (hopefully) reach a consensus genotype call.
- consensus/consensus_called_chr{chr}_{start}-{end}.csv, e.g., consensus_called_chr1_103000029-113000029.csv contains the consensus genotype calls across all datasets. Whenever there is a mismatch, we attempt to find the consensus (majority vote) across all datasets and if found, we replace the genotype call for the strain at a location with the consensus call.
- consensus/no_consensus_chr{chr}_{start}-{end}.csv contains those strains and genomic locations where no consensus was found and they are not compliments. These are removed from the merged dataset and will be imputed in the imputation step.
- consensus/compliments_chr{chr}_{start}-{end}.csv contains the strain and genomic locations where no consensus was found but the mismatched genotype calls are compliments. These are kept in the merged dataset but collected in this file.
- consensus/consensus_summary_chr{chr}_{start}-{end}.csv contains the proportion of agreement with the consensus call and number compared for each dataset. If a dataset has a low degree of agreement with the consensus calls across the other datasets, it should be checked for a dataentry error. These will be summarized across each chromosome and genome (when complete) and should be reviewed.
- summary/summary_table_chr{chr}_{start}-{end}.csv provides an overall all view of how the dataset was merged. It lists the number of strains and bps that were added per dataset through the mergeing process and the final number of bps and strains in the merged dataset.

## IMPUTATION

This pipeline performs imputation using the merged dataset on a given chromosome and region (chr, start, end).

This pipeline performs imputation using the merged dataset on a given chromosome and region (chr, start, end).

For each chromosome region (start, end),

The prepare_for_imputation.R function reads in the merged dataset for the region, calculates phylogenetic distance (if not calculated) and creates the list inputs. Note that phylogenetic distance calculation is computationally expensive so the distance matrix is saved after calculation. The list inputs has length equal to the number of strains in the dataset so each list entry contains information about each strain. If the number of strains > 50, consider distributing the calculation. Otherwise, call get_best_strains() directly.

Now the input files and haploqa python script to run the imputation for each strain is ready.

2. If running in a slurm environment, use generate_array_jobs.R to generate a .sh file that will run imputations in an array job for a chr and region. Once generated, run the .sh file. If not running in a slurm environment, separately run each haplohmm.py in each strain's directory.

3. Once all the haplohmm.py files are complete, run check_haplohmm_output.R to make sure they completed without errors. If there are errors, they will be noted.

4. Once they all run without errors, merge the imputed data with the merged data using get_imputed.R

5. Process (clean)  all results and compute summary statistics with process_imputation_results.R

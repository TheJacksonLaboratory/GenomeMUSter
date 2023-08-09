# GenomeMUSter
GenomeMUSter is available at https://muster.jax.org/ 

GenomeMUSter mouse genetic variation service enables multi-trait, multi-population data integration and analyses

Robyn L. Ball, Molly A. Bogue, Hongping Liang, Anuj Srivastava, David G. Ashbrook, Anna Lamoureux, Matthew W. Gerring, Alexander S. Hatoum, Matthew Kim, Hao He, Jake Emerson, Alexander K. Berger, David O. Walton, Keith Sheppard, Baha El Kassaby, Francisco Castellanos, Govind Kunde-Ramamoorthy, Lu Lu, John Bluis, Sejal Desai, Beth A. Sundberg, Gary Peltz, Zhuoqing Fang, Gary A. Churchill, Robert W. Williams, Arpana Agrawal, Carol J. Bult, Vivek M. Philip, Elissa J. Chesler

Data from sixteen variation datasets (Table 1)  were merged by chromosome and within chromosome in 10 Mb segmented regions. Once the datasets were merged and harmonized, missing genotypes were imputed for each strain in each 10 Mb region.

Code is provided for a segment of the genome, i.e. chr, start, end. For purposes of demonstration, we provide data for a small region.

## MERGE DATASETS

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

For each chromosome or set of chromosomes,
1. Edit prepare_for_imputation_main.R with the chromosome and then call the prepare_for_imputation.R function by submitting prepare_for_imputation.sh. This will result in a number of .sh files, one for each region, in the submit/ directory. Submit these with 'sbatch submit_jobs.sh'

The prepare_for_imputation.R function reads in the merged dataset for the region, calculates phylogenetic distance (if not calculated) and creates the list inputs. Note that phylogenetic distance calculation is computationally expensive so the distance matrix is saved after calculation. The list inputs has length equal to the number of strains in the dataset so each list entry contains information about each strain. If the number of strains > 50, consider distributing the effort of finding the best predictive strains for each target strain in size of 50 blocks (50 strains per block). Otherwise, it calls get_best_strains() directly.

Now the input files and haploqa python script to run the imputation for each strain are ready.

2. When complete, check that all haplohmm.py files are ready for all strains. Edit 'check_prepare_for_imputation_output.R' and set the chromosome.
Run 'sbatch run_check_prepare.sh', look at the log file logs/check_prepare.log to make sure all strains in all regions are ready to impute.
If there are errors, the algorithm will put the missing .sh files in /submit. Run submit_jobs.sh to rerun those jobs and then sbatch run_check_prepare.sh again until there are no more errors.

3. use generate_array_jobs.R to get all imputations to run in an array job for each chr and region. submit all sh files in the submit dir

4. Once all the haplohmm.py files are complete, run check_haplohmm_output.R to make sure they completed without errors. If there are errors, they will be noted on the log file and the .sh files will be moved to the corresponding "../data/out/{chr}/imputation_data/chr{chr}_{start}-{end}/submit" directory. cd to this directory and run sbatch submit_jobs.sh

5. Once they all run without errors, merge the imputed data with the merged data using get_imputation_results_main.R by running run_get_imputation_results_main.sh
Check that the number of files for merged_chr{chr}* and merged_and_imputed_chr{chr}* are the same (no errors)

6. Process all results and compute summary statistics with run_process_imputation_results_{chr}.sh. This script outputs summary statistics across the chromosome for each strain

7. Provide merged_and_imputed_* and imputed_prob_* for ingest to database

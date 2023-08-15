
# Example parameter values
chr <- 1
start <- 103000029
end <- 113000029
indir <- "../data/datasets/"
outdir <- paste0("../data/out/", chr, "/")
strain_idx <- 5

options(stringsAsFactors=FALSE)
source("merge_datasets.R")
source("prepare_for_imputation.R")
source("generate_array_jobs.R")
source("check_haplohmm_output.R")
source("get_imputed.R")
source("process_imputation_results.R")

# strain_idx is the index of the first column in the input files that is a strain. 
#	assumes strain_idx to the last column in the input files are strains
merge_datasets(chr=chr, start=start, end=end, indir=indir, outdir=outdir, strain_idx=5)

# for efficiency, use only 55 of 657 strains in the remainder of the pipeline
df <- fread("../data/out/1/merged_chr1_103000029-113000029.csv", sep=",", data.table=FALSE)
include <- c(1:(strain_idx-1), seq(strain_idx, ncol(df), 12))
df <- df[, include]
imp <- fread("../data/out/1/imputed_prob_chr1_103000029-113000029.csv", sep=",", data.table=FALSE)
imp <- imp[, include]
fwrite(df, "../data/out/1/merged_chr1_103000029-113000029.csv", sep=",", row.names=FALSE)
fwrite(imp, "../data/out/1/imputed_prob_chr1_103000029-113000029.csv", sep=",", row.names=FALSE)

# prepare for imputation by calcualting the phylogenetic distance in the region,
#	identifying the optimal predictor strains for each target strain,
#	set up haplohmm.py input files and directories for each strain
prepare_for_imputation(chr=chr, start=start, end=end, strain_idx=5)
#
# if in a slurm environment, use generate_array_jobs.R to generate a .sh file that
# will run the imputation for all strains. If not in a slurm environment, run
# haplohmm.py for each strain in each region.
#
generate_array_jobs(chr=chr)
# if in a slurm environment, use sbatch to submit the .sh files in the submit/ directory
# once completed, check to see that haplohmm.py completed for each strain. If not, run again.
check_haplohmm_output(chr=chr)
#
# merge imputed genotypes with observed genotypes
# calculate summarty statistics for the region
get_imputed(chr=chr, start=start, end=end, strain_idx=strain_idx, testset=TRUE)
# clean up the data, calculate summary statistics for strains
process_imputation_results(chr=chr, strain_idx=strain_idx, testset=TRUE)



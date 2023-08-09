
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

# strain_idx is the index of the first column in the input files that is a strain. 
#	assumes strain_idx to the last column in the input files are strains
merge_datasets(chr=chr, start=start, end=end, indir=indir, outdir=outdir, strain_idx=5)

prepare_for_imputation(chr=chr, start=start, end=end, strain_idx=5)


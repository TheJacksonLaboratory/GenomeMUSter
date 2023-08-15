# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  for a given strain, creates input files for the HMM (haplohmm.py) from observed data. Copies over haplohmm.py to the strain directory
# INPUT
#       dd:             dataframe (df is default) of observed allele state calls in the region
#	x:		list of inputs for a strain in a region (strain = strain name, d = phylogenetic distance)
#	outdir:		character directory where strain data are stored
#	testset:	logical (default=TRUE) for whether there is a hold out test set
#
# OUTPUT
#                       snps.csv: SNPs that were kept and used in the imputation process
#			snp_removed.csv: indices in dd of the snps that were removed from the imputation step (insuffucuent information)
#			haplotype_ab_codes.csv: input to the HMM (observed data from optimal predictor strains)
#			observation_ab_codes.csv: input to the HMM (observed target strain data)
#			test.csv: set aside test set for this strain and region 
#
#
generate_haploqa_input <- function(dd=df, x=inputs[[i]], outdir=outdir, testset=TRUE) {
  source("df_to_numeric.R")
  source("get_strain_name.R")
  #
  spath <- paste0(outdir, get_strain_name(x$strain), "/")
  subdir <- paste0(outdir, "submit/")
  if (!dir.exists(subdir)) {
    dir.create(subdir)
  }
  #
  print(paste0("Processing ",x$strain))
  x$snp_remove <- as.vector(which(x$alleles[,'major'] == "NA"))
  #dd <- dd[-x$snp_remove, ]
  # create test set for this strain/region
  set.seed(300)
  is_nucleotide <- function(x) {
      x <- as.character(x)
      (!is.na(x) & x!="" & x!="N" & x!="NA")
    }
  if (testset) {
    set.seed(2020)
    known <- which(is_nucleotide(dd[, x$strain]))
    # randomly same 1000 or 10% of the known SNPs (whichever is smaller) and reserve for testing 
    is.test <- sample(known, size = min(1000, round(0.1*length(known))))
    test <- data.frame(chr=rep(dd$chr[1], length(is.test)), bp38=dd$bp38[is.test], test=dd[is.test, x$strain])
    x$test <- test
  }
  #
  if (length(x$snp_remove) > 0)  { 
    x$snps <- dd$bp38[-x$snp_remove]
  } else {
    x$snps <- dd$bp38
  }
  n_snps_removed <- length(x$snp_remove)
  if (length(x$snps) > 0) {
    if (testset) {
      # make all test cases unknown = 0
      dd[is.test, x$strain] <- ""
    }
    if (length(x$snp_remove) > 0)  { 
      x$input_ab_numeric <- df_to_numeric(dd[-x$snp_remove, x$best_strains],
                                                 x$alleles[-x$snp_remove,2:4 ])
    x$obs_numeric <- df_to_numeric(dd[-x$snp_remove, x$strain],
                                             x$alleles[-x$snp_remove,2:4 ])
    } else {
      x$input_ab_numeric <- df_to_numeric(dd[, x$best_strains],
                                                 x$alleles[,2:4 ])
      x$obs_numeric <- df_to_numeric(dd[, x$strain],
                                             x$alleles[,2:4 ])
    }
    #
    write.csv(x$snp_remove, paste0(spath, "snp_remove.csv"), quote=F, row.names = F)
    write.csv(x$snps, paste0(spath, "snps.csv"), quote=F, row.names = F)
    write.csv(x$input_ab_numeric, paste0(spath, "haplotype_ab_codes.csv"), quote=F, row.names = F)
    write.csv(x$obs_numeric, paste0(spath, "observation_ab_codes.csv"), quote=F, row.names = F)
    if (testset) write.csv(x$test, paste0(spath, "test.csv"), quote=F, row.names = F)
    saveRDS(spath, file="inputs.RDS")
    #
    print(paste0("Input files for haploqa saved in: ", spath))
    # copy in code
    dirname <- get_strain_name(x$strain)
    system(command = paste0('cp haplohmm.py ', spath, 'haplohmm.py'))
    shname <- paste0(subdir, 'run_haplohmm_',dirname,'.sh')
    system(command = paste0("sed -e 's/STRAIN/",dirname,"/g' run_haplohmm.sh > ", shname))
    #
  } else {
    print(paste0("Check ", x$strain))
  }
  n_snps_removed
}

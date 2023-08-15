# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  generates a .sh file for slurm enviroment that will run imputation jobs for each strain as an array job
# INPUT
#       chr:              character or numeric
# OUTPUT
#			writes a table in each chromosome region's directory called config.txt with strain names
#			generates a .sh array job file for each region in the chromosome
#
#
generate_array_jobs <- function(chr) {
  options(stringsAsFactors=F)
  cdir <- paste0("../data/out/", chr, "/imputation_data/")
  chr_dirs <- list.files(cdir)
  if (!dir.exists("submit/")) dir.create("submit")
  for (chr_dir in chr_dirs) {
    shname <- paste0("submit/impute_", chr_dir, ".sh")
    system(paste0("cp impute_array_job_main.sh ", shname))
    system(paste0("sed -i 's/CHR_DIR/", chr_dir, "/g' ", shname))
    system(paste0("sed -i 's/CHR/", chr, "/g' ", shname))
    strains <- list.files(paste0(cdir, chr_dir))
    strains <- strains[grep("^[A-Z]|^[1-9]|^i[A-Z]", strains)]
    njobs <- length(strains)
    system(paste0("sed -i 's/NJOBS/", njobs, "/g' ", shname))
    # write strain directories to config.txt
    write.table(paste0(strains, "/"), paste0(cdir, chr_dir, "/config.txt"), col.names=F, row.names=F, quote=F)
  }  
}

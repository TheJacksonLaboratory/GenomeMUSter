# INPUTS
#   chr: chromsome, e.g., 1, 2, ..., X, Y, M
#   start: start bp of the region
#   end: end bp of the region
#   indir: is where the datasets to be merged are stored, organized by dataset name and region, e.g., Sanger_REL2004_chr1_103000029-113000029.csv
#   outdir: is where to write the merged dataset as well as any summary files, mismatches, etc.
#   verbose = TRUE (default) will print the progress of the algorithm. Set to FALSE to silent.
#   strain_idx: the index of the first column that contains strain data
#
# EXAMPLE
#chr <- 1
#start <- 103000029
#end <- 113000029
#indir <- "../data/datasets/"
#outdir <- paste0("../data/out/", chr, "/")
#strain_idx <- 5

merge_datasets <- function(chr, start, end, indir, outdir, strain_idx, verbose=TRUE) {
  # options libraries, and source code (functions used)
  options(stringsAsFactors = FALSE)
  #
  source("match_bps.R")
  source("is_compliment.R")
  source("check_calls.R")
  source("add_strains.R")
  source("add_bps.R")
  source("consensus.R")
  source("get_consensus.R")
  source("get_freq.R")
  source("make_clean.R")
  source("get_agreement.R")
  source("get_ref_alt.R")
  source("check_compliments.R")
  library(data.table)
  #
  if (!dir.exists(outdir)) { dir.create(outdir) }
  #
  mdir <- paste0(outdir, "mismatches/")
  if (!dir.exists(mdir)) { dir.create(mdir) }
  #
  mdir <- paste0(outdir, "consensus/")
  if (!dir.exists(mdir)) { dir.create(mdir) }
  #
  mdir <- paste0(outdir, "summary/")
  if (!dir.exists(mdir)) { dir.create(mdir) }
  #
  dnames <- c("Broad2", "CGD-MDA1", "CGD-MDA2", "CGD-MDA3", "CGD-MDA4", "CGD-MDA5",
              "Perlegen2", "Sanger_REL2004", "Sanger_REL1505","UCLA1", "UNC-GMUGA1", "UNC-MMUGA2", "CC", "BXD", "Stanford1", "B6Eve")
  nfiles <- length(dnames)
  dset <- list()
  dset$name <- dnames
  dset$size <- rep(0, nfiles)
  for (i in 1:length(dnames)) {
    if (verbose) { print(paste0("Extracting SNP data for ", dnames[i])) }
      fn <- paste0(indir, dnames[i], "_chr", chr, "_", start, "-", end, ".csv")
      skip_to_next <- !file.exists(fn)
      if (skip_to_next) { # the file does not exist
        print(paste0("ERROR : No data in ", dnames[i], " for this region. Skipping..."))
        dset$data[[i]] <- NULL
        if (dnames[i] == "B6Eve") {
        # substitute B6 data for B6Eve data if there's no B6Eve data in the region
          print("****Substituting B6 data for B6Eve data in this region\n")
          dset$data[[i]] <- dset$data[[9]][, c(1:6, which(colnames(dset$data[[9]])=="C57BL/6J"))]
          colnames(dset$data[[i]])[7] <- "C57BL/6J Eve"
          cnames <- array(colnames(dset$data[[i]]), dim=c(1,7))
          skip_to_next <- FALSE
        }
      } else { # the file exists
        has_data <- nrow(read.table(fn, sep=",", nrows=2, header=TRUE)) > 0
        if (has_data) { # the dataset is not empty (has rows)
          dset$data[[i]] <- fread(fn, sep=",", data.table=F)
        } else { # the file exists but the dataset is empty
          print(paste0("ERROR : No data in ", dnames[i], " for this region. Skipping..."))
          dset$data[[i]] <- NULL
          skip_to_next <- TRUE
        }
      }
        # get the size of the dataset and the strain names in the dataset
      dset$size[[i]] <- prod(dim(dset$data[[i]]))
      dset$strains[[i]] <- colnames(dset$data[[i]])[strain_idx:ncol(dset$data[[i]])]
  }
  
  nstrains <- length(Reduce(union, dset$strains))
  print(paste0("Loaded ", nstrains, " strains."))

  if ( !all(dset$size == 0) ) {

    # sort so we can load by size to minimize computations
    names(dset$size) <- dset$name
    dset$size <- sort(dset$size, decreasing = TRUE)
    # remove empty datasets
    dset$size <- dset$size[which(dset$size > 0)]
    idx <- match(names(dset$size), dnames)
    # fix list indices to match 
    dset$data <- dset$data[idx]
    dset$name <- dset$name[idx]
    dset$strains <- dset$strains[idx]
    nfiles <- length(dset$data)

    # ready to merge. collect summary information as we merge
    # keep summary in dtable

    dtable <- data.frame(dataset=dset$name, nstrains=rep(NA,nfiles),
                         nlocations=rep(NA,nfiles),
                         strains_added=rep(NA,nfiles),
                         locations_added=rep(NA, nfiles))
    df <- dset$data[[1]]
    df$chr <- gsub("chr", "", df$chr)
    strains <- dset$strains[[1]]
    dtable[1,2:5] <- c(length(strains), nrow(df), length(strains), nrow(df))

    # merge the rest of the datasets

    if (nfiles > 1) {
      for (i in 2:nfiles) {
          nnew <- nnewbp <- 0
          dset$data[[i]]$chr <- gsub("chr", "", dset$data[[i]]$chr)
          # are there strains in these data not included already in the full data
          new.strains <- any(!(dset$strains[[i]] %in% strains))
          # are there bps in these data not included already in the fill data
          new.bps <- any(!(dset$data[[i]]$bp38 %in% df$bp38))
          # add columns for new strains
          if (new.strains) {
            nnew <- length(setdiff(dset$strains[[i]], strains))
            print(paste0('          ',nnew, ' new strains added.'))
            df <- add_strains(addfrom = dset$data[[i]], addto = df,
                            strains = strains, dset_strains=dset$strains[[i]], pname = dset$name[i],
                            checkcalls = TRUE, outdir=outdir, chr=chr, start=start, end=end)
            strains <- union(strains, dset$strains[[i]])
          }
          if (new.bps) {
            nnewbp <- length(setdiff(dset$data[[i]]$bp38, df$bp38))
            print(paste0('          ',nnewbp, ' new bps added.'))
            df <- add_bps(addfrom = dset$data[[i]], addto = df)
          }
          dtable[i,2:5] <- c(length(dset$strains[[i]]), nrow(dset$data[[i]]), nnew, nnewbp)
      }
      # save the summary information
      dtable <- rbind(dtable, c('merged', length(strains), nrow(df), NA, NA))
      fn <- paste0(outdir, 'summary/summary_table_chr', chr, "_", start, "-", end,'.csv')
      write.csv(dtable,fn, quote=F, row.names=F)
     # for any B6 missing, replace with B6Eve. For any B6Eve missing, replace with B6
      if (all(c("C57BL/6J", "C57BL/6J Eve") %in% colnames(df))) {
        b6 <- df[, "C57BL/6J"]
        b6eve <- df[, "C57BL/6J Eve"]
        replace <- which((is.na(b6) | b6 %in% c("", "N")) & b6eve %in% c("A", "T", "C", "G", "H"))
        if (length(replace) > 0) { df[replace, "C57BL/6J"] <- df[replace, "C57BL/6J Eve"] }
        replace <- which((is.na(b6eve) | b6eve %in% c("", "N")) & b6 %in% c("A", "T", "C", "G", "H"))
        if (length(replace) > 0) { df[replace, "C57BL/6J Eve"] <- df[replace, "C57BL/6J"] }
      }
    }

    # clean up the data. sort by genomic location
    df <- df[order(df$bp38), ]
    # clean up mismatches. Get consensus across all datasets for mismatches
    # if it cannot be resolved and it's not a compliment, remove it

    # read in mismatches and get consensus
    pattern <- paste0('mismatches_chr', chr, '_', start, '-', end, '_')
    fns <- list.files(paste0(outdir, "mismatches/"), pattern = pattern)
    mm <- read.csv(paste0(outdir, "mismatches/",fns[1]))
    for (i in 2:length(fns)) {
      mm <- rbind(mm, read.csv(paste0(outdir, "mismatches/",fns[i])))
    }
    result <- get_consensus(mm, dset)
    called <- result$called

    if (!is.null(result$called)) {
      write.csv(result$called, paste0(outdir, "consensus/consensus_called_chr", chr, "_", start, "-", end, ".csv"),
              row.names = FALSE, quote = FALSE)
    }
    if (!is.null(result$mismatches)) {
      write.csv(result$mismatches, paste0(outdir, "consensus/no_consensus_chr", chr, "_", start, "-", end, ".csv"),
              row.names = FALSE, quote = FALSE)
    }
    if (!is.null(result$compliments)) {
      write.csv(result$compliments, paste0(outdir, "consensus/compliments_chr", chr, "_", start, "-", end, ".csv"),
              row.names = FALSE, quote = FALSE)
    }

    if (verbose) {
      print(paste0(length(called$consensus), ' SNP calls will be replaced by the consensus across datasets.'))
      print(paste0('For ', nrow(result$mismatches), ' a consensus could not be reached. These will be removed.'))
    }

    # add consensus calls
    clean <- make_clean(df, called, result$mismatches)
    clean$observed <- gsub("chr", "", clean$observed)
    clean$observed <- apply(clean[, strain_idx:ncol(clean)], 1, get_ref_alt)
    # save it
    print(length(colnames(clean)[strain_idx:ncol(clean)]))
    fn <- paste0(outdir, "merged_chr", chr, "_", start, "-", end, ".csv")
    fwrite(clean, fn, sep=",", row.names = FALSE)
    #
    # add imputed probability matrix, same size as original matrix
    #
    is_nucleotide <- function(x) {
      x <- as.character(x)
      as.numeric(!is.na(x) & x!="" & x!="N")
    }
    imputed_prob <- clean
    imputed_prob[, strain_idx:ncol(imputed_prob)] <- apply(imputed_prob[, strain_idx:ncol(imputed_prob)], 2, is_nucleotide)
    # save imputed prob matrix
    fn <- paste0(outdir, "imputed_prob_chr", chr, "_", start, "-", end, ".csv")
    fwrite(imputed_prob, fn, , sep=",", row.names = FALSE)
    #
    print(paste0("The merged dataset contains ", nrow(clean), " SNPs and ", length(strain_idx:ncol(clean)), " strains."))
  } else {
    if (verbose) {print("No genotype data in this region.")}
  }
}


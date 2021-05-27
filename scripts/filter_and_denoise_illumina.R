#!/usr/local/packages/r-4.0.3/bin/Rscript
options(
    show.error.locations = TRUE,
    show.error.messages = TRUE,
    keep.source = TRUE,
    warn = 1,
    error = function() {
      # cat(attr(last.dump,"error.message"))
      sink(file = stderr())
      dump.frames("dump", TRUE)
      cat('\nTraceback:', file = stderr())
      cat('\n', file = stderr())
      traceback(2) # Print full traceback of function calls with all parameters. The 2 passed to traceback omits the outermost two function calls.
      if (!interactive()) quit(status = 1)
    },
    stringsAsFactors = FALSE
)

require(jsonlite)

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
  dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))
config_file <- file.path(pipelineDir, "config.json")
config <- jsonlite::read_json(
    path = file.path(config_file)
)
.libPaths(config[["r-lib-4.0"]])

library("argparse")
library(dplyr)

parser <- ArgumentParser(description = "Do DADA2 filtering and denoising on Illumina reads")
parser$add_argument(
    "--truncLenL",
    metavar = "BP",
    type = "integer"
)
parser$add_argument(
    "--truncLenR",
    metavar = "BP",
    type = "integer"
)
parser$add_argument(
    "--maxN",
    metavar = "N",
    type = "integer"
)
parser$add_argument(
    "--maxEE",
    metavar = "N",
    type = "character"
)
parser$add_argument(
    "--truncQ",
    metavar = "Q",
    type = "integer"
)
parser$add_argument(
    "--rm.phix",
    metavar = "0|1",
    type = "integer"
)

args <- parser$parse_args()
args$maxEE <- as.numeric(args$maxEE)

library("dada2")
packageVersion("dada2")

run_meta <- function(new_params = list(), checkpoints = list(), samples = list()) {
    info_file <- ".meta.json"

    if (!file.exists(info_file)) {
        run_info <- list()
    } else {
        run_info <- jsonlite::read_json(
            path = info_file
        )
    }

    new_info <- list(params = new_params, checkpoints = checkpoints, samples = samples)
    run_info <- utils::modifyList(run_info, new_info)

    jsonlite::write_json(run_info, info_file, auto_unbox = T)
    invisible(run_info)
}

## perform filtering and trimming
filtpath <- "filtered"

fastqFs <- sort(list.files(pattern="R1_tc.fastq"))
fastqRs <- sort(list.files(pattern="R2_tc.fastq"))
sample.names <- sapply(strsplit(basename(fastqFs), "_"), `[`, 1)
filtF_files<-file.path(filtpath, paste0(sample.names, "_F_filt.fastq.gz"))
filtR_files<-file.path(filtpath, paste0(sample.names, "_R_filt.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
out <- filterAndTrim(
    fwd = fastqFs,
    filt = filtF_files,
    rev = fastqRs,
    filt.rev = filtR_files,
    truncLen = c(args$truncLenL, args$truncLenR),
    maxN = args$maxN,
    maxEE = args$maxEE,
    truncQ = args$truncQ,
    rm.phix = args$rm.phix,
    compress = TRUE,
    multithread = TRUE,
    verbose = TRUE,
    matchIDs = TRUE
)
samples <- lapply(seq_along(fastqFs), function(s) {
    f_fastq <- fastqFs[s]
    # filterAndTrim() returns a matrix with rows named by its fwd argument
    filtF_file <- if (out[f_fastq, "reads.out"] == 0) NULL else filtF_files[[s]]
    filtR_file <- if (out[f_fastq, "reads.out"] == 0) NULL else filtR_files[[s]]
    ans <- list(
        trimmed_fwd = fastqFs[[s]],
        trimmed_rev = fastqRs[[s]],
        filtered_fwd = filtF_file,
        filtered_rev = filtR_file
    )
    return(ans)
}) %>% setNames(sample.names)
run_meta(samples = samples)

## use this to produce the first 5 files to check the quality of the samples.
  ##make post-trimmed quality figures
  ##postscript("filt_forward_reads_quality.eps")
  ##plotQualityProfile(filtFs) ## can take some time if being produced for all samples
  ##dev.off()
  ##if you specify a range it makes only that number of plots, otherwise all plots are produced in grid fashion
  ##postscript("filt_reverse_reads_quality.eps")
  ##plotQualityProfile(filtRs)

  ## Learn errors
  filtFs <- list.files(filtpath, pattern="_F_filt.fastq.gz", full.names = TRUE)
  filtRs <- list.files(filtpath, pattern="_R_filt.fastq.gz", full.names = TRUE)
  sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.namesR
  set.seed(100)
  # Learn forward error rates
  errF <- learnErrors(filtFs, nbases=1e9, multithread=TRUE)
  # Learn reverse error rates
  errR <- learnErrors(filtRs, nbases=1e9, multithread=TRUE)
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
      derepF <- derepFastq(filtFs[[sam]])
      ddF <- dada(derepF, err=errF, multithread=TRUE)
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=errR, multithread=TRUE)
      merger <- mergePairs(ddF, derepF, ddR, derepR)
      mergers[[sam]] <- merger
  }

  rm(derepF); rm(derepR)

  ## Make sequence abundance table 
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, "dada2_abundance_table.rds")

  getN <- function(x) sum(getUniques(x))
  ## track <- cbind(out, rowSums(seqtab))
  v<-rowSums(seqtab)
  v0<-numeric(nrow(out))
  track<-cbind(out, v0)
  rownames(track)<-Map(function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1], rownames(track))
  track[names(v),3]<-v
  colnames(track) <- c("input", "filtered", "merged")
  write.table(
      track,
      "dada2_part1_stats.txt",
      quote = FALSE, 
      append = FALSE, 
      sep = "\t", 
      row.names = TRUE, 
      col.names = TRUE
  )

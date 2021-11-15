#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
  dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))
source(file.path(pipelineDir, "lib", "utils.R"))
require("argparse")

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

library(dplyr)
library("dada2")
packageVersion("dada2")
library(tibble)
library(magrittr)
require(ShortRead)
library(tidyr)

## perform filtering and trimming
filtpath <- "filtered"

fastqFs <- sort(list.files(pattern="R1_tc\\.fastq(?:\\.gz)?"))
fastqRs <- sort(list.files(pattern="R2_tc\\.fastq(?:\\.gz)?"))
sample.names <- sapply(strsplit(basename(fastqFs), "_"), `[`, 1)
filtF_files<-file.path(filtpath, paste0(sample.names, "_F_filt.fastq.gz"))
filtR_files<-file.path(filtpath, paste0(sample.names, "_R_filt.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
nTrimmed <- filterAndTrim(
    fwd = fastqFs,
    filt = filtF_files,
    rev = fastqRs,
    filt.rev = filtR_files,
    maxN = Inf,
    truncQ = 0,
    rm.phix = F,
    compress = F,
    multithread = TRUE,
    verbose = F,
    matchIDs = TRUE
)[, "reads.out", drop = F] %>% set_rownames(sample.names)

nFiltered <- filterAndTrim(
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
)[, "reads.out", drop = F] %>% set_rownames(sample.names)

samples <- lapply(seq_along(fastqFs), function(s) {
    sample <- sample.names[s]
    # filterAndTrim() returns a matrix with rows named by its fwd argument
    filtF_file <- if (nFiltered[sample, "reads.out"] == 0) NULL else filtF_files[[s]]
    filtR_file <- if (nFiltered[sample, "reads.out"] == 0) NULL else filtR_files[[s]]
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
  filtFs <- list.files(filtpath, pattern="_F_filt\\.fastq\\.gz", full.names = TRUE)
  filtRs <- list.files(filtpath, pattern="_R_filt\\.fastq\\.gz", full.names = TRUE)
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

count_reads <- function(file_name) {
    if(file.exists(file_name)) {
        ShortRead::countLines(file_name) / 4
    } else {
        stop("File does not exist")
    }
}

nRaw <- tryCatch({
    # try to use split_library_stats.txt to get stats on raw reads quickly
    read.table("./fwdSplit/split_library_stats.txt", header = T) %>%
        mutate(Input = Reads) %>%
        select(-Reads)
}, error = function(e) {
    # that failed, try counting each sample's raw reads
    tryCatch({
        nRaw <- lapply(gsub(fastqFs, pattern = "_tc\\.fastq", replace = ".fastq"), function(basename) {
            filepath <- list.files(
                path = file.path("fwdSplit", "split_by_sample_out"), 
                pattern = basename,
                full.names = T
            )[1]

            count_reads(filepath)
        }) %>% unlist()

        data.frame(Sample = rownames(nTrimmed), Input = nRaw)
    }, error = function(e) {
        # if that failed, just initialize a data.frame with NA
        data.frame(Sample = rownames(nTrimmed), Input = rep(NA, nrow(nTrimmed)))
    })
}
)

nTrimmed <- nTrimmed %>%
    set_colnames("Trimmed") %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
    
nFiltered <- nFiltered %>%
    set_colnames("Filtered") %>%
    as.data.frame() %>%
    rownames_to_column("Sample")

nMerged <- seqtab %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    mutate(Denoised = rowSums(.[-1])) %>%
    select(Denoised, Sample)

track <- nRaw %>% 
    merge(y = nTrimmed, all.x = T) %>%
    merge(y = nFiltered, all.x = T) %>% 
    merge(y = nMerged, all.x = T) %>%
    dplyr::mutate_at(vars(-Sample), ~tidyr::replace_na(.x, 0)
    ) %>%
    arrange(Sample) %>%
    column_to_rownames("Sample")
  #rownames(track)<-Map(function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1], rownames(track))

  #colnames(track) <- c("Input", "Filtered", "Denoised")
  write.table(
      track,
      "dada2_part1_stats.txt",
      quote = FALSE, 
      append = FALSE, 
      sep = "\t", 
      row.names = TRUE, 
      col.names = TRUE
  )

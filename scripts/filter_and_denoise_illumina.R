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
    action = "store_true"
)
parser$add_argument(
	"--optimize_trim",
	action = "store_true"
)
parser$add_argument(
	"--amplicon_length",
	metavar = "N",
	type = "integer"
)
parser$add_argument(
	"--verbose",
	action = "store_true"
)
parser$add_argument(
	"--debug",
	action = "store_true",
	help = "Use extremely small sample sizes to find errors. Do NOT use this in production execution. Also forces --no-multithread."
)
parser$add_argument(
	"--no-multithread",
	action="store_true",
	help = "Use to speed up small test runs where multithreading would actually slow down execution"
)

args <- parser$parse_args()
args$maxEE <- as.numeric(args$maxEE)
if(args$debug) { args$no_multithread <- T }

library(dplyr)
library("dada2")
packageVersion("dada2")
library(tibble)
library(magrittr)
require(ShortRead)
library(tidyr)
library(Biostrings)

count_reads <- function(file_name) {
    if(file.exists(file_name)) {
        ShortRead::countLines(file_name) / 4
    } else {
        stop("File does not exist")
    }
}

## perform filtering and trimming
filtpath <- "filtered"

fastqFs <- sort(list.files(pattern="R1_tc\\.fastq(?:\\.gz)?"))
fastqRs <- sort(list.files(pattern="R2_tc\\.fastq(?:\\.gz)?"))
sample.names <- sapply(strsplit(basename(fastqFs), "_"), `[`, 1)
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filt<-file.path(filtpath, paste0(sample.names, "_F_tc.fastq.gz"))
filt.rev<-file.path(filtpath, paste0(sample.names, "_R_tc.fastq.gz"))
names(filt) <- names(filt.rev) <- sample.names

learn_errors_and_denoise <- function(filt, filt.rev, fast = F, verbose = T, debug = F, multithread = T) {
	# Learn error rates
	nBases = if(debug) 10^2 else if(fast) 10^6 else 10^9
	if(!verbose) { s <- file(tempfile()); sink(file=s)}
	errF <- learnErrors(
		filt, nbases=nBases, multithread=multithread, randomize = T, verbose = verbose,
		MAX_CONSIST= if(debug) 2 else 10
	)
	errR <- learnErrors(
		filt.rev, nbases=nBases, multithread=multithread, randomize = T, verbose = verbose,
		MAX_CONSIST= if(debug) 2 else 10
	)
	if(!verbose) { sink(); close(s) }
	
	# Sample inference and merger of paired-end reads
	derepF <- derepFastq(filt)
	ddF <- dada(derepF, err=errF, multithread=multithread, verbose = F)
	derepR <- derepFastq(filt.rev)
	ddR <- dada(derepR, err=errR, multithread=multithread, verbose = F)
	merger <- mergePairs(ddF, derepF, ddR, derepR)
	rm(derepF); rm(derepR)

	## Make sequence abundance table 
	seqtab <- makeSequenceTable(merger)
	return(seqtab)
}

eval_trim_L <- function(truncLen_R, fwd, rev, filt, filt.rev, args) {
	# truncLen_R becomes a fixed value in the function returned
    f <- function(truncLen_L) {
        reads.out <- suppressMessages(filterAndTrim(
            fwd = fwd,
            filt = filt,
            rev = rev,
            filt.rev = filt.rev,
            truncLen = c(round(truncLen_L), round(truncLen_R)),
            maxN = args$maxN,
            maxEE = args$maxEE,
            truncQ = args$truncQ,
            rm.phix = args$rm.phix,
            compress = F,
            multithread = !(args$no_multithread),
            verbose = FALSE,
            matchIDs = TRUE
        )[, "reads.out", drop = T])

		filt <- filt[reads.out > 0]
		filt.rev <- filt.rev[reads.out > 0]

		seqtab <- learn_errors_and_denoise(
			filt = filt, filt.rev = filt.rev, fast = T, verbose = args$verbose,
			debug = args$debug, multithread = !(args$no_multithread)
		)
		if(args$verbose) {
			message(paste("Params:", truncLen_L, truncLen_R))
			message(paste("Value:", sum(seqtab)))
		}
		return(sum(seqtab))
    }
	return(f)
}
eval_trim_R <- function(truncLen_L, fwd, rev, filt, filt.rev, args) {
	# truncLen_R becomes a fixed value in the function returned
    f <- function(truncLen_R) {
        reads.out <- suppressMessages(filterAndTrim(
            fwd = fwd,
            filt = filt,
            rev = rev,
            filt.rev = filt.rev,
            truncLen = c(round(truncLen_L), round(truncLen_R)),
            maxN = args$maxN,
            maxEE = args$maxEE,
            truncQ = args$truncQ,
            rm.phix = args$rm.phix,
            compress = F,
            multithread = !(args$no_multithread),
            verbose = FALSE,
            matchIDs = TRUE
        )[, "reads.out", drop = T])

		filt <- filt[reads.out > 0]
		filt.rev <- filt.rev[reads.out > 0]

		seqtab <- learn_errors_and_denoise(
			filt = filt, filt.rev = filt.rev, fast = T, verbose = args$verbose,
			debug = args$debug, multithread = (!(args$no_multithread))
		)
		if(args$verbose) {
			message(paste("Params:", truncLen_L, truncLen_R))
			message(paste("Value:", sum(seqtab)))
		}
        return(sum(seqtab))
    }
	return(f)
}

optimize_trim_len <- function(fwd, rev, filt, filt.rev, args) {
	min_total_len <- args$amplicon_len + 12

	max_len_L <- max(nchar(readDNAStringSet(fwd[1], format = "fastq")))
	max_len_R <- max(nchar(readDNAStringSet(rev[1], format = "fastq")))

	first_time <- T
	maxit <- if (! args$debug) 10 else 2
	curit <- 1
	trunc_lens <- list(L = max_len_L, R = max_len_R)
	old_vals <- c(trunc_lens$L, trunc_lens$R, 0)
	seqtab <- NULL
	while(curit < maxit) {
		r_func <- eval_trim_R(
			trunc_lens$L, fwd = fastqFs,
            filt = filt,
            rev = fastqRs,
            filt.rev = filt.rev, 
			args = args)
		lower <- max(
			60, min_total_len - max_len_L, min_total_len - trunc_lens$L
		)
		if(args$verbose){
			message(paste("Optimizing truncLen_R. Min/max:", lower, max_len_R))
		}
		o <- optim(
			fn = r_func, 
			par = trunc_lens$R, 
			lower = lower,
			upper = max_len_R,
			method="Brent", 
			control = list(
				reltol = if(! args$debug) 3 else (max_len_R - lower) * 0.6, 
				fnscale = -1
				)
			)
		if(!first_time) {
			if(round(o$par) == old_vals[2] || o$value == old_vals[3]) {
				break()
			}
		}
		trunc_lens$R <- round(o$par)
		if(args$verbose) {
			message(paste("Current F length:", trunc_lens$L))
			message(paste("Current R length:", trunc_lens$R))
			message(paste("Num merged:", o$value))
		}
		old_vals <- c(trunc_lens$L, trunc_lens$R, o$value)

		f_func <- eval_trim_L(trunc_lens$R, fwd = fastqFs,
            filt = filt,
            rev = fastqRs,
            filt.rev = filt.rev, 
			args = args			)
		lower <- max(60, min_total_len - max_len_R, min_total_len - trunc_lens$R)
		if(args$verbose){
			message(paste("Optimizing truncLen_L. Min/max:", lower, max_len_L))
		}
		o <- optim(
			fn = f_func, 
			par = trunc_lens$L, 
			lower = lower,
			upper = max_len_L,
			method="Brent", 
			control = list(
				reltol = if(! args$debug) 3 else (max_len_L - lower) * 0.6, 
				fnscale = -1
				)
			)
		if(!first_time) {
			if(round(o$par) == old_vals[1]) {
				break()
			}
		}
		trunc_lens$L <- round(o$par)
		if(args$verbose) {
			message(paste("Current F length:", trunc_lens$L))
			message(paste("Current R length:", trunc_lens$R))
			message(paste("Num merged:", o$value))
		}
		old_vals <- c(trunc_lens$L, trunc_lens$R, o$value)

		first_time <- F
		curit <- curit + 1
	}

	print(paste("Optimal F length:", trunc_lens$L))
	print(paste("Optimal R length:", trunc_lens$R))

	return(list(L = trunc_lens$L, R = trunc_lens$R))
}

if(args$optimize_trim) {
	trunc_lens <- optimize_trim_len(fwd = fastqFs,
            filt = filt,
            rev = fastqRs,
            filt.rev = filt.rev, args = args)
	# after calling this function, the trimmed files are all present...but we
	# still don't have separate trimming and filtering stats
	args$truncLenL <- trunc_lens$L
	args$truncLenR <- trunc_lens$R
}

nTrimmed <- filterAndTrim(
    fwd = fastqFs,
    filt = filt,
    rev = fastqRs,
    filt.rev = filt.rev,
    maxN = Inf,
    truncQ = 0,
    rm.phix = F,
    compress = F,
    multithread = !(args$no_multithread),
    verbose = F,
    matchIDs = TRUE
)[, "reads.out", drop = F] %>% set_rownames(sample.names)

nFiltered <- filterAndTrim(
    fwd = fastqFs,
    filt = filt,
    rev = fastqRs,
    filt.rev = filt.rev,
    truncLen = c(args$truncLenL, args$truncLenR),
    maxN = args$maxN,
    maxEE = args$maxEE,
    truncQ = args$truncQ,
    rm.phix = args$rm.phix,
    compress = TRUE,
    multithread = !(args$no_multithread),
    verbose = TRUE,
    matchIDs = TRUE
)[, "reads.out", drop = F] %>% set_rownames(sample.names)

samples <- lapply(seq_along(fastqFs), function(s) {
    sample <- sample.names[s]
    # filterAndTrim() returns a matrix with rows named by its fwd argument
    filtF_file <- if (nFiltered[sample, "reads.out"] == 0) NULL else filt[[s]]
    filtR_file <- if (nFiltered[sample, "reads.out"] == 0) NULL else filt.rev[[s]]
    ans <- list(
        trimmed_fwd = fastqFs[[s]],
        trimmed_rev = fastqRs[[s]],
        filtered_fwd = filtF_file,
        filtered_rev = filtR_file
    )
    return(ans)
}) %>% setNames(sample.names)
run_meta(samples = samples)

filt <- filt[nFiltered[, "reads.out"] > 0]
filt.rev <- filt.rev[nFiltered[, "reads.out"] > 0]
if(args$optimize_trim) {
	# use this to produce the first 5 files to check the quality of the samples.
	# make post-trimmed quality figures
	postscript("filt_forward_reads_quality.eps")
	plotQualityProfile(filt) ## can take some time if being produced for all samples
	dev.off()
	# if you specify a range it makes only that number of plots, otherwise all plots are produced in grid fashion
	postscript("filt_reverse_reads_quality.eps")
	plotQualityProfile(filt.rev)
	dev.off()
}

seqtab <- learn_errors_and_denoise(filt = filt, filt.rev = filt.rev, debug = args$debug, multithread = !(args$no_multithread))
saveRDS(seqtab, "dada2_abundance_table.rds")

# Collect statistics
nRaw <- tryCatch({
    # try to use split_library_stats.txt to get stats on raw reads quickly
    read.table("./fwdSplit/split_library_stats.txt", header = T) %>%
        mutate(Input = Reads) %>%
        select(-Reads)
}, error = function(e) {
    # that failed, try counting each sample's raw reads
    tryCatch({
        nRaw <- lapply(gsub(fastqFs, pattern = "_tc\\.fastq(?:\\.gz)", replace = ".fastq"), function(basename) {
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

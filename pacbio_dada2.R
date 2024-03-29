#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
  dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
source(file.path(pipelineDir, "lib", "utils.R"))
require("argparse")

parser <- ArgumentParser(
    description = "MSL 16S pipeline for preprocessing and denoising PacBio runs",
    epilog =
        paste(
            "This script is analogous to the MSL 16S pipeline for Illumina",
            "runs. Reads will be trimmed of adapters, quality-trimmed,",
            "quality-filtered, and denoised using DADA2 functions."
        )
)
parser$add_argument(
    "input",
    metavar = "INPUT_DIRECTORY",
    type = "character",
    help = "Directory containing demuxed files (.ccs.fastq.gz)"
)
parser$add_argument(
    "--pattern",
    metavar = "GLOB",
    type = "character",
    help = "Regex pattern that will match all DEMUXED .ccs.fastq.gz files. The first capture group will be extracted as sample name. Remember to escape regex special characters. For correct bash syntax, enclose the pattern in double quotes."
)
parser$add_argument(
    "--wd",
    metavar = "PATH",
    type = "character",
    help = "Working directory. The directory's base name will be taken as the run name.",
    required = T
)
parser$add_argument(
    "--nodelete",
    "--no-delete",
    "--no_delete",
    action = "store_true",
    help = "Don't delete intermediate sequence files after denoising."
)
parser$add_argument(
    "--forward_primer",
	default = "AGRGTTYGATYMTGGCTCAG",
	type = "character",
	action = "store",
    help = "Forward primer to remove. Sequences without both primers are discarded."
)
parser$add_argument(
    "--reverse_primer",
	default = "RGYTACCTTGTTACGACTT",
	type = "character",
	action = "store",
    help = "Reverse primer to remove. Sequences without both primers are discarded."
)
parser$add_argument(
    "--max_mismatch",
	type = "integer",
	action = "store",
	default = 2,
    help = "The number of mismatches to tolerate when matching reads to primer sequences."
)
parser$add_argument(
    "--indels",
	action = "store_true",
	default = T,
    help = "Allow insertions or deletions of bases when matching adapters"
)
# the DADA2 tutorial default for max_length is 1600, but 4000 is the max
# seq length in our SILVA db
parser$add_argument(
    "--max_length",
	action="store", default=4000,
    type = "integer",
	help = "Remove reads with length longer than max_length. max_length is enforced before quality trimming and truncation."
)
parser$add_argument(
    "--min_length",
	action="store", default=1000,
    type = "integer",
	help = "Remove reads with length shorter than min_length. min_length is enforced after quality trimming and truncation."
)
parser$add_argument(
    "--rm.phix",
    action = "store_true"
)
parser$add_argument(
	"--memory",
	type = "integer",
	default = 64,
	help = "(GB) Limits the number of threads dada2 functions are allowed to start, based on an internal model of memory used per thread. The default is sufficient for PacBio Sequel II runs."
)
parser$add_argument(
	"--multithread",
	default = parallel::detectCores(),
	action="store",
	type = "integer",
	help = "Manually set number of threads. This parameter will be overridden if parallel::detectCores() finds fewer cores."
)
# the command line options are responsible for telling this script what
# resources are /available/. This script then decides how much of those 
# resources to use.

args <- parser$parse_args()
if (any(is.null(args))) {
    stop("Some args missing!")
}
args$multithread <- max(min(parallel::detectCores(), args$multithread), 1)

setwd(args$wd)
run <- basename(getwd())
sink(
    file = file.path(getwd(), paste0(run, "_16S_pipeline_log.txt")), split = T, append = T
)
run_dir <- getwd()

require("dada2")
packageVersion("dada2") # 1.12.1 

require("ShortRead")
library(dplyr)
library(magrittr)
library(tidyr)

# share all output files if we're working in the global run directory
if(grepl(run_dir, pattern = paste0("^", config[["run_storage_path"]]))) {
    Sys.umask(mode=0002)
}

inPath <- args$input

cache_checksums <- function(files, name) {
    checksums <- lapply(files, function(f) {
        gsub("^\\s+|\\s+$", "", system(paste("sha512sum", f), intern = T))
    }) %>% setNames(files)

    checksums_list <- list()
    checksums_list[[name]] <- checksums

    run_meta(checkpoints = checksums_list)
    return()
}
run_meta(new_params = list(platform = "PACBIO"))

# glob_pattern <- gsub("()", "(.*?)", args$pattern, fixed = T)
files <- list.files(inPath, full.names = T)
fastqs <- sort(files[grepl(args$pattern, files, perl = T)])
print(fastqs)
fastqs_basename <- sapply(fastqs, basename)

rex <- regexpr(args$pattern, fastqs_basename, perl=TRUE)
sample.names <- paste(
	run, 
	substr(
		fastqs_basename, 
		attr(rex, 'capture.start')[,"sample"], 
		attr(rex, 'capture.start')[,"sample"] + attr(rex, 'capture.length')[,"sample"] - 1
	), 
	sep = "."
	)

input_size <- sum(file.size(fastqs)) # bytes

# copy to a local directory regardless of original location, and pre-pend the
# run name while we're at it
demuxed_path <- file.path("demultiplexed")
dir.create(demuxed_path, showWarnings = T)
demuxed_filepaths <- file.path(
    demuxed_path, paste(sample.names, "fastq.gz", sep = ".")
    )
targets <- !fastqs %in% list.files(demuxed_filepaths, full.names = T)
file.copy(
    fastqs[targets], demuxed_filepaths[targets],
    overwrite = T, copy.mode = TRUE, copy.date = FALSE
)
cache_checksums(demuxed_filepaths, "samples")

names(demuxed_filepaths) <- sample.names
samples <- data.frame(CCS = demuxed_filepaths)

trim_primers <- function(ins, outs) {
    if (length(ins) == 0) {
        cat("No samples to trim primers from\n")
        return()
    }

	arg_list <- list(
		fn = ins,
		fout = outs,
		primer.fwd = args$forward_primer,
		primer.rev = args$reverse_primer,
		max.mismatch = args$max_mismatch, 
		allow.indels = args$indels,
		orient = T, verbose = T
	)
	print("Executing dada2::removePrimers() with args:")
	print(arg_list)

    # Remove primers. Discard reads without primers.
    prim.stats <- do.call(
						dada2::removePrimers,
						arg_list
					) %>%
					set_rownames(names(ins))
	
	# class(prim.stats) [1] "matrix" "array", so use `[` slicing operator
    if(all(prim.stats[, "reads.out"] == 0)) {
		warning("No reads passed the Removing Primers step  (Did you select the right primers?)")
	}
	print("")
    
	# get stats AND delete the output file if removePrimers reports no reads passing
    samples <- lapply(seq_along(ins), function(s) {
        fastq <- names(ins)[s]
		if (prim.stats[fastq, "reads.out"] == 0 && file.exists(outs[[s]])) {
			file.remove(outs[s])
		}
        # removePrimers returns a matrix with rows named by its fwd argument
        trimmed_file <- if (prim.stats[fastq, "reads.out"] == 0) NULL else outs[[s]]
        ans <- list(
            raw = ins[[s]],
            trimmed = trimmed_file
        )
        
        return(ans)
    }) %>% setNames(names(ins))
    run_meta(samples = samples)

    return(prim.stats)
}

trim_and_filter <- function(ins, outs) {
    if (length(ins) == 0) {
        cat("No samples to QC")
        return()
    }
    # cat("Quality-trimming and filtering: \n")
    # if (!is.null(names(ins))) {
    #     cat(paste(names(ins), "\n"))
    # } else {
    #     cat(paste(ins, "\n"))
    # }

    targets <- if (length(ins) > 5) {
        sample(x = 1:length(ins), size = 5, replace = FALSE)
    } else {
        seq_along(ins)
    }
    png("01_tagcleaned_reads_quality.png", width = 1020, height = 600)
    print(suppressMessages(dada2::plotQualityProfile(ins[targets])))
    dev.off()

    png("02_tagcleaned_reads_quality_aggregate.png", width = 1020, height = 600)
    print(dada2::plotQualityProfile(ins, aggregate = TRUE))
    dev.off()

    quality_trimmed <- paste0(outs, ".qtrim.fastq.gz")

	input_size <- sum(file.size(ins))
	mt <- max(
		min(
			args$multithread, 
			floor((args$memory-12.7) / (3.61 * input_size/1E9 -2.6))
		), 
		1
	)

	arg_list <- list(
		fwd = ins,
        filt = quality_trimmed,
		truncQ = 2,
        minQ = 0,
        maxN = Inf,
        maxEE = Inf,
		minLen = 0,
        rm.phix = args$rm.phix,
        compress = TRUE,
        multithread = mt,
        verbose = TRUE
	)

	print("Executing dada2::filterAndTrim() with args:")
	print(arg_list)

	gc()

    ## just count read lengths after trimming N's
	do.call(dada2::filterAndTrim, arg_list)

    lens.fn <- lapply(quality_trimmed[file.exists(quality_trimmed)], function(fn) {
        nchar(dada2::getSequences(fn))
    })
    lens <- do.call(c, lens.fn)

    png("03_quality_trimmed_sample_lens.png",
        width = 4743,
        height = 400
    )
    hist(lens, lwd = 1, breaks = 100)
    abline(v = args$min_length, col = "#E69F00")
    abline(v = args$max_length, col = "#E69F00")
    abline(h = 0)
    dev.off()

	arg_list$filt <- outs
	arg_list$maxLen <- args$max_length
	arg_list$minQ <- 3
	arg_list$maxEE <- 2
	arg_list$minLen <- args$min_length
	cat(paste("Filtering with max pre-trim length:", args$max_length, "\n"))
	cat(paste("Filtering with min post-trim length:", args$min_length, "\n"))

    # redo filter and trim, this time with length filtering
    filter_in_out_counts <-	do.call(dada2::filterAndTrim, arg_list)
    
    samples <- lapply(seq_along(ins), function(s) {
        # filterAndTrim() returns a matrix with rows named by its fwd argument
        if(! names(ins)[s] %in% rownames(filter_in_out_counts) || filter_in_out_counts[names(ins)[s], "reads.out"] == 0) {
            filt_file <- NULL
        } else {
            filt_file <- outs[[s]]
        }
        ans <- list(
            trimmed = ins[[s]],
            filtered = filt_file
        )
        
        return(ans)
    }) %>% setNames(names(ins))
    run_meta(samples = samples)
    
    
    file.remove(quality_trimmed)
    return(filter_in_out_counts)
}

denoise <- function(ins) {
    cat("Denoising all samples... \n")

    set.seed(100)

    # dereplicate
	arg_list <- list(
		fls = ins,
		qualityType = "FastqQuality",
		n = 5e+05,
		verbose = T
	)
	print("Executing dada2::derepFastq() with args:")
	print(arg_list)
    drp <- do.call(dada2::derepFastq, arg_list)
	
    # Learn error rates
	arg_list <- list(
		BAND_SIZE = 32,
		multithread = args$multithread,
		errorEstimationFunction = dada2:::PacBioErrfun,
        randomize = T,
		verbose = T
	)
	print("Executing dada2::learnErrors() with args:")
	print(arg_list)
	err <- do.call(dada2::learnErrors, c(list(fls = drp), arg_list))
	
    png("04_error_profiles.png", height = 1000, width = 1000)
    dada2::plotErrors(err)
    dev.off()

    # Sample inference
	arg_list <- list(
	)
    print("Executing dada2::dada() with args:")
	print(arg_list)
	dd <- do.call(dada2::dada, c(list(derep = drp, err = err), arg_list))
	
    rm(drp)

	print("")

    ## Make sequence abundance table
    return(dada2::makeSequenceTable(dd))
}

collectStats <- function(
    in_files = list(), 
    primer_trimmed_files = list(), 
    filtered_files = list(),
    denoised_table = data.frame()
) {
    cat("Collecting stats...\n")

    count_reads <- function(file_name) {
        if(file.exists(file_name)) {
            ShortRead::countLines(file_name) / 4
        } else {
            return(0)
        }
    }
    in_counts <- sapply(in_files, count_reads)
    primer_trimmed_counts <- sapply(primer_trimmed_files, count_reads)
    filtered_counts <- sapply(filtered_files, count_reads)
    
    denoised_counts <- data.frame(denoised = rowSums(denoised_table), sample_name = rownames(denoised_table))

    stats <- cbind.data.frame(in_counts, primer_trimmed_counts, filtered_counts)
    stats$sample_name <- names(in_files)
    stats <- stats %>% 
        merge(denoised_counts, all.x = T) %>%
        select(-sample_name) %>%
        replace(is.na(.), 0) %>%
        set_rownames(names(in_files))

    return(stats)
}

tagcleanedpath <- file.path("tagcleaned")
dir.create(tagcleanedpath, showWarnings = TRUE)
tcs <- file.path(tagcleanedpath, paste0(sample.names, "_trimmed.fastq.gz"))
names(tcs) <- sample.names
samples <- cbind(samples, primer_trimmed = tcs)
targets <- !tcs %in% list.files(tagcleanedpath, full.names = T) & file.exists(samples$CCS)
# don't process repeats

primer_trim_output <- trim_primers(
    ins = samples$CCS[targets] %>% setNames(rownames(samples)[targets]),
    outs = samples$primer_trimmed[targets]
)

filtpath <- file.path("filtered")
dir.create(filtpath, showWarnings = TRUE)
filtereds <- (paste0(sample.names, "_filt.fastq.gz"))
filtereds_files <- file.path(filtpath, filtereds)
names(filtereds_files) <- names(tcs)
samples <- cbind(samples, filtered = filtereds_files)
targets <- !filtereds_files %in% list.files(filtpath, full.names = T) & file.exists(samples$primer_trimmed)

trim_and_filter_output <- trim_and_filter(
    ins = samples$primer_trimmed[targets] %>% 
        setNames(rownames(samples)[targets]),
    outs = samples$filtered[targets]
)

seq_tab_output <- "dada2_abundance_table.rds"
seq_tab <- denoise(filtereds_files[file.exists(filtereds_files)])
saveRDS(seq_tab, seq_tab_output)

seq_tab <- readRDS(seq_tab_output)
stats <- collectStats(
    in_files = samples$CCS %>% setNames(rownames(samples)),
    primer_trimmed_files = samples$primer_trimmed,
    filtered_files = samples$filtered,
    denoised_table = seq_tab
)
colnames(stats) <- c("Input", "Trimmed", "Filtered", "Denoised")
write.table(
    stats, "dada2_part1_stats.txt",
    quote = FALSE, append = FALSE, sep ="\t", row.names = TRUE, col.names = TRUE
)

if(! args$nodelete){
    unlink("tagcleaned/*_trimmed.fastq.gz", expand = T)
    if(length(list.files("tagcleaned/")) == 0)  {
        unlink("tagcleaned", recursive = T)
    } else {
        warning("Couldn't remove directory ./tagcleaned/ because it contained unrecognized files.")
    }

    unlink("filtered/*_filt.fastq.gz", expand = T)
    if(length(list.files("filtered/")) == 0)  {
        unlink("filtered", recursive = T)
    } else {
        warning("Couldn't remove directory ./filtered/ because it contained unrecognized files.")
    }
}

remove_chimeras <- function(args) {
    bim2 <- dada2::isBimeraDenovo(
        seqtab,
        minFoldParentOverAbundance = 3.5, multithread = args$multithread
    )
    table(bim2) # FALSE: 1794; TRUE: 1290
    sum(seqtab[, bim2]) / sum(seqtab) # 0.04069179
    # Lots of ASVs labeled as bimeras, but only make up 4% of reads
}
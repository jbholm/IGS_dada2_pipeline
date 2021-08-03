#!/usr/bin/env Rscript
options(
    show.error.locations = TRUE,
    show.error.messages = TRUE,
    keep.source = TRUE,
    warn = 1,
    error = function() {
        # cat(attr(last.dump,"error.message"))
        sink(file = stderr())
        dump.frames("dump", TRUE)
        cat("\nTraceback:", file = stderr())
        cat("\n", file = stderr())
        traceback(2) # Print full traceback of function calls with all parameters. The 2 passed to traceback omits the outermost two function calls.
        if (!interactive()) quit(status = 1)
    },
    stringsAsFactors = FALSE
)

require(jsonlite)

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
  dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
config_file <- file.path(pipelineDir, "config.json")
config <- jsonlite::read_json(
    path = file.path(config_file)
)

.libPaths(config[["r-lib"]])
require("dada2")
packageVersion("dada2") # 1.12.1 
require("argparse")
require("ShortRead")
library(dplyr)
library(magrittr)
library(tidyr)

parser <- ArgumentParser(description = "Assign taxa")
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
    help = "Regex pattern that will match all DEMUXED .ccs.fastq.gz files. Use () in place of sample name. Remember to escape regex special characters. For correct bash syntax, enclose the pattern in double quotes."
)
parser$add_argument(
    "--wd",
    metavar = "PATH",
    type = "character",
    help = "Working directory. The directory's base name will be taken as the run name."
)
args <- parser$parse_args()
if (any(is.null(args))) {
    stop("Some args missing!")
}
setwd(args$wd)
run <- basename(getwd())
sink(
    file = file.path(getwd(), paste0(run, "_16S_pipeline_log.txt")), split = T, append = T
)

run_dir <- getwd()
# share all output files if we're working in the global run directory
if(grepl(run_dir, pattern = paste0("^", config[["run_storage_path"]]))) {
    Sys.umask(mode=0002)
}

inPath <- args$input

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

glob_pattern <- gsub("()", "(.*)", args$pattern, fixed = T)
fastqs <- sort(list.files(inPath, pattern = glob_pattern, full.names = T)) # B01\..+\.css\.fastq\.gz
print(fastqs)
sample.names <- sapply(fastqs, function(filename) {
    paste(run, sub(glob_pattern, "\\1", basename(filename), perl = T), sep = ".")
})

# copy to a local directory regardless of original location, and pre-pend the
# run name while we're at it
demuxed_path <- file.path("demuxed")
dir.create(demuxed_path, showWarnings = T)
demuxed_filepaths <- file.path(
    demuxed_path, paste(sample.names, "fastq.gz", sep = ".")
    )
file.copy(
    fastqs, demuxed_filepaths,
    overwrite = T, copy.mode = TRUE, copy.date = FALSE
)
cache_checksums(demuxed_filepaths, "samples")


names(demuxed_filepaths) <- sample.names
samples <- data.frame(CCS = demuxed_filepaths)

trim_primers <- function(ins, outs) {
    F27 <- "AGRGTTYGATYMTGGCTCAG" # Pacbio primers given in DADA2 tutorials
    R1492 <- "RGYTACCTTGTTACGACTT"

    if (length(ins) == 0) {
        cat("No samples to trim primers from")
        return()
    }

    cat("Trimming primers from: \n")
    if (!is.null(names(ins))) {
        cat(paste(names(ins), "\n"))
    } else {
        cat(paste(ins, "\n"))
    }

    # Remove primers. Discard reads without primers. Write these counts somewhere?
    prim.stats <- t(mapply(function(ins, outs) {
        stats <- tryCatch(
            {
                dada2::removePrimers(
                    ins,
                    outs,
                    primer.fwd = F27,
                    primer.rev = dada2:::rc(R1492),
                    orient = TRUE
                )
            },
            error = function(e) {
                fq <- readFastq(ins)
                inseqs <- length(fq)
                ans <- data.frame(reads.in = inseqs, reads.out = 0)
                return(ans)
            }
        )
        return(stats)
    }, ins = ins, outs = outs)) %>%
        set_rownames(names(ins)) %>%
        set_colnames(c("reads.in", "reads.out"))
    print(prim.stats)
    
    samples <- lapply(seq_along(ins), function(s) {
        fastq <- names(ins)[s]
        # filterAndTrim() returns a matrix with rows named by its fwd argument
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
    cat("Quality-trimming and filtering: \n")
    if (!is.null(names(ins))) {
        cat(paste(names(ins), "\n"))
    } else {
        cat(paste(ins, "\n"))
    }

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
    ## perform quality trimming
    dada2::filterAndTrim(
        fwd = ins,
        filt = quality_trimmed,
        truncQ = 3,
        minQ = 0,
        maxN = Inf,
        maxEE = Inf,
        rm.phix = FALSE, # DADA2 PacBio tutorials did not remove PhiX
        compress = TRUE,
        multithread = T,
        verbose = TRUE
    )

    lens.fn <- lapply(quality_trimmed[file.exists(quality_trimmed)], function(fn) {
        nchar(dada2::getSequences(fn))
    })
    lens <- do.call(c, lens.fn)

    # min_length <- 0
    # max_length <- max(lens)
    # min_length_chosen <- max_length_chosen <- F
    # prev_thresh <- 0
    # auto_thresh <- -1
    # current_lens <- lens
    # bad_thresh <- F
    # while (!min_length_chosen && !max_length_chosen && !bad_thresh) {

    #     auto_thresh <- autothresholdr::auto_thresh(current_lens, "Renyi")[1] # 1589bp
    #     if (auto_thresh > min(current_lens) && auto_thresh < 1200) {
    #     # SILVA tells us there are 16S genes as short as 1200. Be open to novel
    #     # sequences
    #     min_length_chosen <- T
    #     min_length <- auto_thresh
    #     bad_thresh <- F
    #     } else if (auto_thresh > 1600 && auto_thresh < max(current_lens)) {
    #     # 1600 is B Callahan's default. SILVA tells us there are even longer ones
    #     # out there
    #     max_length_chosen <- T
    #     max_length <- auto_thresh
    #     bad_thresh <- F
    #     } else {
    #     bad_thresh <- T
    #     }
    #     current_lens <- current_lens[current_lens >= min_length & current_lens <= max_length]
    # }
    # if (!min_length_chosen) {
    #     min_length <- 1000
    # } else {
    #     warning(paste0("Minimum read length decreased from 1000 to ", min_length))
    # }
    # if (!max_length_chosen) {
    #     max_length <- 1600 # defaults in Pacbio tutorials
    # } else {
    #     warning(paste0("Maximum read length increased from 1600 to ", auto_thresh))
    # }

    min_length <- 1000
    max_length <- 4000 # the DADA2 tutorial default is 1600, but 4000 is the max
    # seq length in our SILVA db
    png("03_quality_trimmed_sample_lens.png",
        width = 4743,
        height = 400
    )
    hist(lens, lwd = 1, breaks = 100)
    abline(v = min_length, col = "#E69F00")
    abline(v = max_length, col = "#E69F00")
    abline(h = 0)
    dev.off()

    # redo filter and trim, this time with length filtering
    filter_in_out_counts <- dada2::filterAndTrim(
        fwd = ins,
        filt = outs,
        minQ = 3,
        maxN = 0,
        maxEE = 1,
        truncQ = 0,
        minLen = min_length,
        maxLen = max_length,
        rm.phix = FALSE,
        compress = TRUE,
        multithread = T,
        verbose = TRUE
    )
    
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
    drp <- dada2::derepFastq(ins, verbose = TRUE, qualityType = "FastqQuality")
    # Learn error rates
    err <- dada2::learnErrors(drp,
        BAND_SIZE = 32, multithread = TRUE,
        errorEstimationFunction = dada2:::PacBioErrfun,
        randomize = TRUE,
        nbases = 1e9
    ) # seconds

    png("04_error_profiles.png", height = 1000, width = 1000)
    dada2::plotErrors(err)
    dev.off()

    # Sample inference
    dd <- dada2::dada(drp, err = err, multithread = TRUE)
    rm(drp)

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
targets <- !tcs %in% list.files(tagcleanedpath, full.names = T)
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
colnames(stats) <- c("Input", "Adapter-trimmed", "Filtered", "Denoised")
write.table(
    stats, "dada2_part1_stats.txt",
    quote = FALSE, append = FALSE, sep = , row.names = TRUE, col.names = TRUE
)

remove_chimeras <- function(args) {
    bim2 <- dada2::isBimeraDenovo(
        seqtab,
        minFoldParentOverAbundance = 3.5, multithread = TRUE
    )
    table(bim2) # FALSE: 1794; TRUE: 1290
    sum(seqtab[, bim2]) / sum(seqtab) # 0.04069179
    # Lots of ASVs labeled as bimeras, but only make up 4% of reads
}
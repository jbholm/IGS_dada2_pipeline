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
    dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))
config_file <- file.path(pipelineDir, "config.json")
config <- jsonlite::read_json(
    path = file.path(config_file)
)
.libPaths(config[["r-lib-4.0"]])

require("argparse")

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <- dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))

parser <- ArgumentParser(description = "Remove chimeras and assign taxa")
parser$add_argument("--seq",
    metavar = "SEQUENCING_MACHINE", type = "character",
    help = "ILLUMINA or PACBIO"
)
parser$add_argument("runs",
    metavar = "RUN_ID", type = "character",
    nargs = "+",
    help = "Name of run directory(ies)"
)
args <- parser$parse_args()
if (is.null(args$seq)) {
    stop("Please provide --seq. Execute dada2_part2.R --help for help.")
}
if (any(!dir.exists(args$runs))) {
    stop(do.call(base::paste, args = list(c(
        "Directories do not exist:",
        args$runs[!dir.exists(args$runs)]
    ), sep = " ")))
}

suppressMessages(
    suppressWarnings({
        library("Biostrings")
        require("dada2")
        library(tidyr)
        library(dplyr)
        library(magrittr)
        library(readr)
    })
)
runs <- args$runs %>% setNames(lapply(args$runs, basename))

project_meta <- function(runs = list(), samples = list()) {
    filepath <- ".meta.json"
    if (!file.exists(filepath)) {
        proj_info <- list()
    } else {
        proj_info <- jsonlite::read_json(
            path = filepath
        )
    }

    if (length(runs) > 0 || length(samples) > 0) {
        proj_info$runs[names(runs)] <- runs
        proj_info$samples <- samples

        jsonlite::write_json(proj_info, filepath, auto_unbox = T)
    }
    return(proj_info)
}

# combine maps
# map_filepath <- (function(runs) {
#     # get each map, looking in the path specified by the run's metadata file
#     maps <- lapply(names(runs), function(run_name) {
#         run_info <- jsonlite::read_json(file.path(runs[run_name], ".meta.json"))
#         if ("map" %in% names(run_info[["checkpoints"]])) {
#             map_filepath <- file.path(
#                 runs[run_name], names(run_info[["checkpoints"]][["map"]])[1]
#             )
#             tryCatch(
#                 {
#                     suppressMessages(
#                         map <- read_tsv(map_filepath, quote = "", na = "", trim_ws = F)
#                     )
#                 },
#                 error = function(e) {
#                     stop(paste("Unable to find mapping file; supposed to be at", map_filepath))
#                 }
#             )

#             if (length(runs) > 1) {
#                 # add a Run column to the map
#                 map <- map %>%
#                     mutate(Run = run_name) %>%
#                     relocate(Run)
#             }
#             return(map)
#         } else {
#             return()
#         }
#     })
#     bind_rows(maps) %>% write_tsv("map.txt")
#     return("map.txt")
# })(runs)

## INPUT
## list all of the files matching the pattern
counts_and_stats <- (function(runs) {
    tables <- lapply(names(runs), function(run_name) {
        count_table <- readRDS(
            list.files(
                runs[run_name],
                pattern = "dada2_abundance_table.rds", full.names = TRUE
            )[[1]]
        )
        # if (length(runs) > 1) {
        #     rownames(count_table) <- paste(
        #         run_name, rownames(count_table),
        #         sep = "_"
        #     )
        # }

        return(count_table)
    })
    stats <- lapply(names(runs), function(run_name) {
        stat_table <- read.csv(
            list.files(
                runs[run_name], 
                pattern = "dada2_part1_stats.txt", 
                full.names = TRUE
                )[[1]],
            sep = "", stringsAsFactors = FALSE
        )
        # if (length(runs) > 1) {
        #     old_sample_names <- rownames(stat_table)
        #     rownames(stat_table) <- paste(
        #         run_name, rownames(stat_table),
        #         sep = "_"
        #     )
        #     run_metadata <- project_meta()$runs[[run_name]]
        #     names(run_metadata$samples) <- rownames(stat_table)[match(names(run_metadata$samples), old_sample_names)]
        #     project_meta(runs = list(run_metadata) %>% set_names(run_name))
        # }

        return(stat_table)
    })
    # basically, this code is here solely to accommodate when two runs have the
    # same name, or a run with two samples of the same name (if that's even possible)
    # old_stats_sample_names <- do.call(c, lapply(stats, rownames))
    # old_count_table_sample_names <- do.call(c, lapply(tables, rownames))

    # new_sample_names <- make.names(old_stats_sample_names, unique = T)
    stats <- bind_rows(stats)
    # stats <- stats %>% set_rownames(new_sample_names)

    samples <- do.call(c, lapply(names(runs), function(run_name) {
        run_metadata <- project_meta()$runs[[run_name]]
        return(run_metadata$samples)
    }))
    # names(samples) <- rownames(stats)[match(names(samples), old_stats_sample_names)]
    project_meta(samples = samples)

    # Combine count tables
    unqs <- unique(c(sapply(tables, colnames), recursive = TRUE))
    n <- sum(unlist(lapply(X = tables, FUN = nrow)))
    st <- matrix(0L, nrow = n, ncol = length(unqs)) %>% set_colnames(unqs)
    row <- 1
    for (table in tables) {
        st[row:(row + nrow(table) - 1), colnames(table)] <- table
        row <- row + nrow(table)
    }

    # Make count table sample names unique (and matching up with stat sample names)
    # stats_i <- 1
    # rownames(st) <- lapply(
    #     seq_along(old_count_table_sample_names),
    #     function(count_table_i) {
    #         while (old_count_table_sample_names[count_table_i] != old_stats_sample_names[stats_i]) {
    #             stats_i <<- stats_i + 1
    #         }
    #         return(new_sample_names[stats_i])
    #     }
    # ) %>%
    #     unlist()

    # sort ASVs by total frequency
    st <- st[, order(colSums(st), decreasing = TRUE)]
    return(list(counts = st, stats = stats))
})(runs)

## st.all<-mergeSequenceTables(runs)
# Remove chimeras
mfpoa <- if (args$seq == "ILLUMINA") 2 else 3.5
seqtab <- dada2::removeBimeraDenovo(counts_and_stats$counts, method = "consensus", minFoldParentOverAbundance = mfpoa, multithread = TRUE)
# Write to disk
saveRDS(seqtab, "all_runs_dada2_abundance_table.rds") # CHANGE ME to where you want sequence table saved
write.csv(seqtab, "all_runs_dada2_abundance_table.csv", quote = FALSE)
outputs <- c("all_runs_dada2_abundance_table.rds", "all_runs_dada2_abundance_table.csv")
# Create a reference FASTA containing all the ASVs from colnames(seqtab)
# (each sequence's FASTA header is itself)
fasta <- "all_runs_dada2_ASV.fasta"
writeXStringSet(DNAStringSet(colnames(seqtab)), fasta, format = "fasta")
outputs <- append(outputs, fasta)

nonchimerics <- rowSums(seqtab)[rownames(counts_and_stats$stats)] %>%
    replace_na(0) %>%
    setNames(rownames(counts_and_stats$stats))
project_stats <- cbind(Sample = rownames(counts_and_stats$stats), counts_and_stats$stats, Nonchimeric = nonchimerics)
write.table(project_stats, "DADA2_stats.txt",
    quote = FALSE, append = FALSE,
    sep = "\t", row.names = F, col.names = TRUE
)
outputs <- append(outputs, c("DADA2_stats.txt", map_filepath))

cat(outputs)
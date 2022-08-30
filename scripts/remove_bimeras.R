#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
    dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))
source(file.path(pipelineDir, "lib", "utils.R"))
require("argparse")


parser <- ArgumentParser(description = "Remove chimeras and assign taxa")
parser$add_argument("--seq",
    metavar = "SEQUENCING_MACHINE", type = "character",
    help = "ILLUMINA or PACBIO"
)
parser$add_argument("--map",
    nargs = "?",
    metavar = "PROJECT_MAP", type = "character",
    help = "Tab-delimited file with two columns: RUN.PLATEPOSITION and sampleID"
)
parser$add_argument("runs",
    metavar = "RUN_ID", type = "character",
    nargs = "+",
    help = "Name of run directory(ies)"
)
parser$add_argument("--verbose", action = "store_true", help = "Enables all STDERR from R functions.")
parser$add_argument("--log", type = "character", help = "Redirect messages and warnings to a file instead of STDERR")

args <- parser$parse_args()
if (is.null(args$seq)) {
    stop("Please provide --seq. Execute scripts/remove_bimeras.R --help for help.")
}
if (any(!dir.exists(args$runs))) {
    stop(do.call(base::paste, args = list(c(
        "Directories do not exist:",
        args$runs[!dir.exists(args$runs)]
    ), sep = " ")))
}
if (!is.null(args$log)) {
    con <- file(args$log, open = "a")
    sink(con, append = TRUE, type = "message")
}

suppress_if_not_verbose({
    library("Biostrings")
    require("dada2")
    library(tidyr)
    library(dplyr)
    library(magrittr)
    library(readr)
})

runs <- args$runs %>% setNames(lapply(args$runs, basename))

## INPUT
## list all of the files matching the pattern
counts_and_stats <- (function(runs) {
    message("Reading in count tables")
    tables <- lapply(names(runs), function(run_name) {
        count_table <- readRDS(
            list.files(
                runs[run_name],
                pattern = "dada2_abundance_table.rds", full.names = TRUE
            )[[1]]
        )

        return(count_table)
    })

    message("Reading in stats tables")
    stats <- lapply(names(runs), function(run_name) {
        stat_table <- read.csv(
            list.files(
                runs[run_name],
                pattern = "dada2_part1_stats.txt",
                full.names = TRUE
            )[[1]],
            sep = "", stringsAsFactors = FALSE
        )

        return(stat_table)
    })

    message("Combining stats")
    stats <- bind_rows(stats)

    message("Recording metadata from runs")
    samples <- do.call(c, lapply(names(runs), function(run_name) {
        run_metadata <- project_meta()$runs[[run_name]]
        return(run_metadata$samples)
    }))
    project_meta(samples = samples)

    # Combine count tables
    message("Combining count tables")
    unqs <- unique(c(sapply(tables, colnames), recursive = TRUE))
    n <- sum(unlist(lapply(X = tables, FUN = nrow)))
    st <- matrix(0L, nrow = n, ncol = length(unqs)) %>%
        set_colnames(unqs) %>%
        set_rownames(c(sapply(tables, rownames), recursive = TRUE))
    row <- 1
    for (table in tables) {
        st[row:(row + nrow(table) - 1), colnames(table)] <- table
        row <- row + nrow(table)
    }

    # sort ASVs by total frequency
    message("Sorting ASVs by total frequency")
    st <- st[, order(colSums(st), decreasing = TRUE)]
    return(list(counts = st, stats = stats))
})(runs)

seqtab <- counts_and_stats$counts
stats <- counts_and_stats$stats

# If given a project map, apply new sample names and remove samples not in map
# trim_ws=T trims whitespace to avoid user errors in map
if (!is.null(args$map)) {
    message("Reading in project map")

    map <- suppress_if_not_verbose({
        read_tsv(
            args$map,
            quote = "",
            na = "",
            trim_ws = T,
            col_types = cols_only(
                RUN.PLATEPOSITION = col_character(), sampleID = col_character()
            )
        )
    })
    plates <- c("A", "B", "C", "D")
    control_positions <- suppress_if_not_verbose(
        read_csv(file.path(pipelineDir, "share", "controls_platepositions.csv"), quote = "", na = "", trim_ws = T)
    )

    message("Looking for UDI plates")

    for (plate in plates) {
        matches <- regexpr(text = map$RUN.PLATEPOSITION, pattern = paste0("UDI", plate, "\\.[A-H]\\.[1-9]{1,2}$"), perl = T)
        if (any(matches > -1)) {
            message(paste("Found instances of UDI plate", plate))

            runs_containing_plate <- map %>%
                mutate(Matches = matches) %>%
                filter(Matches > -1) %>%
                mutate(Run = substr(RUN.PLATEPOSITION, start = 1, stop = Matches - 2)) %>%
                pull(Run) %>%
                unique()

            for (run in runs_containing_plate) {
                message(paste("Checking if controls need to be added for plate", plate, "in run", run))

                controls <- control_positions %>%
                    filter(Plate == plate) %>%
                    mutate(RUN.PLATEPOSITION = paste(run, PLATEPOSITION, sep = ".")) %>%
                    select(RUN.PLATEPOSITION, sampleID)
                controls_to_add <- controls[!controls$RUN.PLATEPOSITION %in% map$RUN.PLATEPOSITION, ]
                if (nrow(controls_to_add) > 0) {
                    message("Adding controls")
                    map <- map %>%
                        rows_insert(controls_to_add, by = "RUN.PLATEPOSITION")
                } else {
                    message("Controls already in project map.")
                }
            }
        }
    }

    message("Writing final project map")
    write_tsv(map, args$map)

    merge_with_map <- function(df) {
        return(
            df %>%
                mutate(RUN.PLATEPOSITION = rownames(df)) %>%
                merge(map) %>%
                mutate(sampleID = make.names(sampleID, unique = T)) %>%
                set_rownames(.$sampleID)
        )
    }
    # write this here for the time-being
    message("Writing temporary all-runs abundance table")
    write.csv(seqtab, "all_runs_dada2_abundance_table.csv", quote = FALSE)

    message("Using map to subset sample data from abundance table")
    seqtab.df <- as.data.frame(seqtab)
    seqtab <- seqtab.df %>%
        merge_with_map() %>%
        select(-c(RUN.PLATEPOSITION, sampleID)) %>%
        as.matrix() %>%
        extract(, apply(., MAR = 2, function(col) sum(col) > 0))

    message("Using map to subset sample data from stats table")
    stats <- counts_and_stats$stats %>%
        merge_with_map() %>%
        select(-c(RUN.PLATEPOSITION, sampleID))
    if (nrow(stats) == 0) {
        stop("Provided map did not match any samples in runs")
    }
    if (nrow(seqtab) == 0) {
        stop("None of the samples in the project map passed denoising.")
    }
}

# Remove chimeras
message("Removing chimeras")
mfpoa <- if (args$seq == "ILLUMINA") 1.5 else 3.5
seqtab <- dada2::removeBimeraDenovo(seqtab, method = "consensus", minFoldParentOverAbundance = mfpoa, multithread = TRUE)

# Write to disk
message("Writing final abunadnce table, stats, ASVs")
saveRDS(seqtab, "all_runs_dada2_abundance_table.rds") # CHANGE ME to where you want sequence table saved
write.csv(seqtab, "all_runs_dada2_abundance_table.csv", quote = FALSE)
outputs <- c("all_runs_dada2_abundance_table.rds", "all_runs_dada2_abundance_table.csv")
# Create a reference FASTA containing all the ASVs from colnames(seqtab)
# (each sequence's FASTA header is itself)
fasta <- "all_runs_dada2_ASV.fasta"
writeXStringSet(DNAStringSet(colnames(seqtab)), fasta, format = "fasta")
outputs <- append(outputs, fasta)

nonchimerics <- rowSums(seqtab)[rownames(stats)] %>%
    replace_na(0) %>%
    setNames(rownames(stats))

project_stats <- cbind(Sample = rownames(stats), stats, Nonchimeric = nonchimerics)

write.table(project_stats, "DADA2_stats.txt",
    quote = FALSE, append = FALSE,
    sep = "\t", row.names = F, col.names = TRUE
)
outputs <- append(outputs, "DADA2_stats.txt")

cat(outputs)
#!/usr/local/packages/r-3.6.0/bin/Rscript
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
.libPaths("/home/jolim/share/R/x86_64-pc-linux-gnu-library/3.6")
require("argparse")
library("Biostrings")

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <- dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))

parser <- ArgumentParser(description='Remove chimeras and assign taxa')
parser$add_argument('--seq', metavar='SEQUENCING_MACHINE', type="character",
                  help='ILLUMINA or PACBIO')
parser$add_argument("runs",
  metavar = "RUN_ID", type = "character",
  nargs = "+",
  help = "Name of run directory(ies)"
)
args <- parser$parse_args()
if (is.null(args$seq)) {
  stop("Please provide --seq. Execute dada2_part2.R --help for help.")
}
if(any(!dir.exists(args$runs))) {
  stop(do.call(base::paste, args = list(c(
    "Directories do not exist:",
    args$runs[!dir.exists(args$runs)]
  ), sep = " ")))
}

require("dada2")
library(tidyr)
library(dplyr)
library(magrittr)

runs <- args$runs %>% setNames(lapply(args$runs, basename))
## INPUT
## list all of the files matching the pattern 
counts_and_stats <- (function(runs) {
    tables<-lapply(runs, function(run) {
        readRDS(list.files(run, pattern = "dada2_abundance_table.rds", full.names = TRUE)[[1]])
    })
    stats <- lapply(runs, function(run) {
        read.csv(list.files(run, pattern = "dada2_part1_stats.txt", full.names = TRUE)[[1]], sep = "", stringsAsFactors = FALSE)
    })
    old_stats_sample_names <- do.call(c, lapply(stats, rownames))
    old_count_table_sample_names <- do.call(c, lapply(tables, rownames))
    
    new_sample_names <- make.names(old_stats_sample_names, unique = T)
    stats <- bind_rows(stats) %>% `row.names<-`(new_sample_names)
    # in future, use bind_rows to create stats using id column, then concatenate
    # run name and rowname to create unique rownames. These can then be used 
    # in place of new_sample_names to name the rows of st

    # Combine count tables
    unqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
    n<-sum(unlist(lapply(X=tables, FUN = nrow)))
    st <- matrix(0L, nrow = n, ncol = length(unqs)) %>% set_colnames(unqs)
    row <- 1
    for (table in tables) {
        st[row:(row + nrow(table) - 1), colnames(table)] <- table
        row <- row + nrow(table)
    }

    # Make count table sample names unique (and matching up with stat sample names)
    stats_i <- 1
    rownames(st) <- lapply(
        seq_along(old_count_table_sample_names), 
        function(count_table_i) {
            while (old_count_table_sample_names[count_table_i] != old_stats_sample_names[stats_i]) {
                stats_i <<- stats_i + 1
            }
            return(new_sample_names[stats_i])
        }
    ) %>%
    unlist()

    # sort ASVs by total frequency
    st <- st[, order(colSums(st), decreasing = TRUE)]
    return(list(counts = st, stats = stats))
})(runs)

##st.all<-mergeSequenceTables(runs)
# Remove chimeras
mfpoa <- if(args$seq == "ILLUMINA") 2 else 3.5
seqtab <- dada2::removeBimeraDenovo(counts_and_stats$counts, method="consensus", minFoldParentOverAbundance=mfpoa, multithread=TRUE)
# Write to disk
saveRDS(seqtab, "all_runs_dada2_abundance_table.rds") # CHANGE ME to where you want sequence table saved
write.csv(seqtab, "all_runs_dada2_abundance_table.csv", quote=FALSE)
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
outputs <- append(outputs, "DADA2_stats.txt")

cat(outputs)

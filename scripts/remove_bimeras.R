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
require("dada2")
library(tidyr)
path <- getwd()

runs <- args$runs %>% setNames(lapply(args$runs, basename))
## INPUT
## list all of the files matching the pattern 
counts_and_stats <- (function(runs) {
  tables<-sapply(runs, function(run) {
    subdir <- file.path(path, run)
    list.files(subdir, pattern="dada2_abundance_table.rds", full.names=TRUE)[[1]]
  })
  stats<-sapply(runs, function(run) {
    subdir <- file.path(path, run)
    list.files(subdir, pattern="dada2_part1_stats.txt", full.names=TRUE)[[1]]
  })


  ## get the run names using splitstring on the tables where - exists
  runstats <- runs <- vector("list", length(tables)) %>% setNames(names(runs))
  for (run in names(runs)) {
    runs[[run]] <- readRDS(tables[run])
    runstats[[run]] <- read.delim(stats[run])
  }

  stats <- do.call(rbind, lapply(stats, function(f) {
    read.csv(f, sep="", stringsAsFactors=FALSE)
  }))

  unqs <- unique(c(sapply(runs, colnames), recursive=TRUE))
  n<-sum(unlist(lapply(X=runs, FUN = nrow)))
  st <- matrix(0L, nrow=n, ncol=length(unqs))
  rownames(st) <- c(sapply(runs, rownames), recursive=TRUE)
  colnames(st) <- unqs
  for(sti in runs) {
    st[rownames(sti), colnames(sti)] <- sti
  }
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
fc = file("all_runs_dada2_ASV.fasta")
fltp = character()
for( i in 1:ncol(seqtab))
{
  fltp <- append(fltp, paste0(">", colnames(seqtab)[i]))
  fltp <- append(fltp, colnames(seqtab)[i])
}
writeLines(fltp, fc)
rm(fltp)
close(fc)
outputs <- append(outputs, "all_runs_dada2_ASV.fasta")

nonchimerics <- rowSums(seqtab)[rownames(counts_and_stats$stats)] %>%
  replace_na(0) %>%
  setNames(rownames(counts_and_stats$stats))
project_stats <- cbind(counts_and_stats$stats, nonchimeric = nonchimerics)
write.table(project_stats, "DADA2_stats.txt",
  quote = FALSE, append = FALSE,
  sep = "\t", row.names = TRUE, col.names = TRUE
)
outputs <- append(outputs, "DADA2_stats.txt")

cat(outputs)

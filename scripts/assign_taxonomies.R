#!/usr/local/packages/r-3.6.0/bin/Rscript
options(
  show.error.locations = TRUE,
  show.error.messages = TRUE,
  error = quote(quit(status = 1))
)
.libPaths("/home/jolim/share/R/x86_64-pc-linux-gnu-library/3.6")
require("argparse")


initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
  dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))

parser <- ArgumentParser(description = 'Assign taxa')
parser$add_argument(
  'input',
  metavar = "FASTA",
  type = "character",
  help = paste0("A FASTA with ASV names in the headers"
                ),
                nargs = 1
  )
parser$add_argument(
  '--tax',
  metavar = 'TAXONOMY',
  type = "character",
  help = "SILVA128, SILVA132, SILVA138, HOMD, or UNITE",
  action = 'append',
  choices = c("SILVA128", "SILVA132", "SILVA138", "HOMD", "UNITE")
)
parser$add_argument(
  '--minBoot',
  metavar = 'BOOTSTRAP_CONFIDENCE',
  type = "integer",
  default = 80
)
args <- parser$parse_args()
# if (is.null(args$input)) {
#   stop("Please provide an ASV table stored in an RDS.\n")
# }
# if (length(args$tax) == 0) {
#   stop("No taxonomy given.\n")
# }
require("dada2")
require("Biostrings")
require("rjson")
path <- getwd()

config_file <- file.path(pipelineDir, "config.json")
config_text <- readChar(config_file, file.info(config_file)$size)
config <- rjson::fromJSON(config_text, simplify = TRUE)

# seqtab <- readRDS(args$input)
# Assign taxonomy (requires colnames of seqtab to be ASV sequences)
taxonomies <-
  c(
    SILVA128 = "silva_nr_v128_train_set.fa.gz",
    SILVA132 = "silva_nr_v132_train_set.fa.gz",
    SILVA138 = "silva_nr99_v138_wSpecies_train_set.fa.gz",
    HOMD = "HOMD_v15.1_DADA2_taxonomy_final.txt",
    UNITE = "sh_general_release_dynamic_01.12.2017.fasta"
  )

fasta <- args$input[[1]]
seqs <- dada2::getSequences(fasta)
dummy <- lapply(args$tax, function(taxonomy) {
  db <- taxonomies[taxonomy]
  db <- file.path(config$ref_16s_dir, db)

  if (taxonomy == "SILVA128") {
    taxLevels = c("Domain",
                  "Phylum",
                  "Class",
                  "Order",
                  "Family",
                  "Genus",
                  "Species") # SILVA doesn't have species, but this is trivial to DADA2
  } else if (taxonomy %in% c("HOMD", "UNITE", "SILVA132", "SILVA138")) {
    taxLevels = c("Kingdom",
                  "Phylum",
                  "Class",
                  "Order",
                  "Family",
                  "Genus",
                  "Species")
  } else {
    stop(
      paste(
        "Runtime error: Unknown taxonomy: ",
        taxonomy,
        ". Please specify taxon levels in code."
      )
    )
  }

  tryCatch({
    assignments <- dada2::assignTaxonomy(seqs,
                                db,
                                multithread = TRUE,
                                taxLevels = taxLevels, outputBootstraps = TRUE, minBoot=args$minBoot)
  }, error = function(e) {
    print("Error in dada2::assignTaxonomy:")
    stop(e)
  })
  assignments$tax <- gsub('^[a-zA-Z]__', '', assignments$tax, perl = TRUE)
  # assignments <- apply(assignments, MAR = 1, function(taxon) {
  #   # remove any taxon level hints (k__, p__, c__, o__, g__, s__)
  #   return(gsub('^[a-zA-Z]__', '', taxon, perl = TRUE))
  # })
  write.csv(assignments$tax,
            paste(taxonomy, "classification.csv", sep = "."),
            quote = FALSE)
  cat(paste(taxonomy, "classification.csv", sep = "."))
})

#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
    dirname(dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)])))
source(file.path(pipelineDir, "lib", "utils.R"))
require("argparse")
if(length(grep("--verbose", initial.options, fixed = T)) > 0) {
    require("dada2")
} else {
	suppressMessages(suppressWarnings(require("dada2")))
}

parser <- ArgumentParser(description = "Assign taxa")
parser$add_argument(
    "input",
    metavar = "FASTA",
    type = "character",
    help = paste0("A FASTA with ASV names in the headers"),
    nargs = 1
)
parser$add_argument(
    "--tax",
    metavar = "TAXONOMY",
    type = "character",
    help = do.call(
        paste, c(as.list(names(config[["taxonomy_dbs"]])), sep = ", ")
    ),
    action = "append",
    choices = names(config[["taxonomy_dbs"]])
)
parser$add_argument(
    "--minBoot",
    metavar = "BOOTSTRAP_CONFIDENCE",
    type = "integer",
    default = formals(dada2::assignTaxonomy)$minBoot
)
parser$add_argument(
    "--verbose",
    required = F, action = "store_true", help = "Enables all STDERR from R functions."
)
parser$add_argument("--log", type = "character", help = "Redirect messages and warnings to a file instead of STDERR")

args <- parser$parse_args()
if (!is.null(args$log)) {
    con <- file(args$log, open = "a")
    sink(con, append = TRUE, type = "message")
}

# if (is.null(args$input)) {
#   stop("Please provide an ASV table stored in an RDS.\n")
# }
# if (length(args$tax) == 0) {
#   stop("No taxonomy given.\n")
# }
suppress_if_not_verbose({
    require("Biostrings")
})
path <- getwd()

# Assign taxonomy (requires colnames of seqtab to be ASV sequences)
taxonomies <- config[["taxonomy_dbs"]]

message("Reading in sequences")
fasta <- args$input[[1]]
seqs <- dada2::getSequences(fasta)

stdout <- lapply(args$tax, function(taxonomy) {
    db <- taxonomies[taxonomy]
    db <- file.path(config$ref_16s_dir, db)
    message(paste(taxonomy, "database:", db))

    if (taxonomy == "SILVA128") {
        taxLevels <- c(
            "Domain",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species"
        ) # SILVA doesn't have species, but this is trivial to DADA2
    } else if (taxonomy %in% c("HOMD", "UNITE", "SILVA132", "SILVA138", "SILVA138forPB")) {
        taxLevels <- c(
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species"
        )
    } else {
        stop(
            paste(
                "Runtime error: Unknown taxonomy:",
                taxonomy,
                ". Please specify taxon levels in code."
            )
        )
    }

    message("Assigning taxonomy")
    tryCatch(
        {
            invisible(
                capture.output(
                    assignments <- dada2::assignTaxonomy(
                        seqs,
                        db,
                        multithread = TRUE,
                        taxLevels = taxLevels,
                        outputBootstraps = TRUE,
                        minBoot = args$minBoot
                    )
                )
            )
        },
        error = function(e) {
            print("Error in dada2::assignTaxonomy:")
            stop(e)
        }
    )
    assignments$tax <- gsub("^[a-zA-Z]__", "", assignments$tax, perl = TRUE)
    # assignments <- apply(assignments, MAR = 1, function(taxon) {
    #   # remove any taxon level hints (k__, p__, c__, o__, g__, s__)
    #   return(gsub('^[a-zA-Z]__', '', taxon, perl = TRUE))
    # })
    rownames(assignments$tax) <- names(seqs)
    write.csv(assignments$tax,
        paste(taxonomy, "classification.csv", sep = "."),
        quote = FALSE
    )
    return(paste0(taxonomy, ".", "classification.csv"))
})
cat(do.call(paste, c(stdout, sep = "\n")))
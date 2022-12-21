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
    default = 80 # DADA2 default is 50, which is low; if doing post-processing, you can lower this default so more of the output is retained
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
    } else if (taxonomy %in% names(config[["taxonomy_dbs"]])) {
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
			sink(file = stderr(), type = "output")
            taxtab <- dada2::assignTaxonomy(
                seqs,
                db,
                multithread = TRUE,
                taxLevels = taxLevels,
                outputBootstraps = TRUE,
                minBoot = args$minBoot,
                verbose = T
            )

            binom <- dada2::assignSpecies(
                seqs, 
                file.path(config$ref_16s_dir, config$species_inhouse_db),
                verbose = T
            )
			sink()
        },
        error = function(e) {
            print("Error in dada2::assignTaxonomy or dada2::assignSpecies:")
            stop(e)
        }
    )
    taxtab$tax <- gsub("^[a-zA-Z]__", "", taxtab$tax, perl = TRUE)
    
    # for each seq where the genuses don't match, delete assignments above genus
    # and override genus with the genus from assignSpecies
    if ("Genus" %in% colnames(taxtab$tax)) {
        gcol <- which(colnames(taxtab$tax) == "Genus")
    } else {
        gcol <- ncol(taxtab)
    }
    # if genuses don't match, and there IS an assignment to species level by
    # exact match, then trust the exact match
    gen.match <- mapply(dada2:::matchGenera, taxtab$tax[, gcol], binom[, 1])
    binom_override <- !gen.match & !is.na(binom[, 2]) 
    taxtab$tax[binom_override, 1:(gcol - 1)] <- NA_character_
    taxtab$tax[binom_override, gcol] <- binom[binom_override, 1]
    
    # Where there was species annotation from an exact match, use it
    species <- binom[, 2]
    if(ncol(taxtab$tax) > gcol ) {
        species[is.na(species)] <- taxtab$tax[is.na(species), gcol + 1]
        taxtab$tax[, gcol + 1] <- species
    } else {
        taxtab$tax <- cbind(taxtab$tax, species)
    }
    rownames(taxtab$tax) <- names(seqs)
    write.csv(taxtab$tax,
        paste(taxonomy, "classification.csv", sep = "."),
        quote = FALSE
    )
    return(paste0(taxonomy, ".", "classification.csv"))
})
cat(do.call(paste, c(stdout, sep = "\n")))
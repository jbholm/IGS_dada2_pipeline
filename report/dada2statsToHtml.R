#!/usr/bin/env Rscript

require(jsonlite)
initial.options <- commandArgs(trailingOnly = FALSE)
pipelineDir <-
    dirname(
        dirname(
            sub("--file=", "", initial.options[grep("--file=", initial.options)])
        )
    )
config_file <- file.path(pipelineDir, "config.json")
config <- jsonlite::read_json(
    path = file.path(config_file)
)
.libPaths(config[["r-lib"]])

require("argparse", quietly = T)

parser <- ArgumentParser(description = "")
parser$add_argument(
    "--file", "-f",
    metavar = "DADA2_STATS",
    type = "character",
    help = "The DADA2 stats file, tab-delimited"
)
parser$add_argument(
    "--group_by", "-g",
    metavar = "COLUMN_INDEX",
    type = "integer",
    help = "Which column index to group by (optional)"
)
opt <- parser$parse_args()

require("kableExtra")
require("stringi")
require("jsonlite")

if(is.null(opt$file)) {
    stop("Statistics file must be provided using -f or --file")
}

dada2 <- read.delim(opt$file, row.names=1)

capitalize <- function(y) {
  vapply(y, function(x) {
    c <- strsplit(x, " ")[[1]]
    paste(toupper(substring(c, 1,1)), substring(c, 2),
          sep="", collapse=" ")
  }, FUN.VALUE = "")
}
colnames(dada2) = capitalize(colnames(dada2))

require(kableExtra)
library(magrittr)

fileConn<-file("dada2stats.json")
#output <- dada2 %>% kableExtra::kable() %>% kableExtra::kable_styling() %>% kableExtra::scroll_box(width = "100%", height = "100%")
#output <- stringi::stri_replace_all(output, "<div class=\"kable\" style", regex = "^<div style")
output <- jsonlite::toJSON(dada2)
writeLines(output, fileConn)
close(fileConn)
#
# discrete <- is.character(mapping[, opt$group_by, drop = TRUE]) || is.factor(mapping[, opt$group_by, drop = TRUE])
# if(discrete) {
#     nGroups <- length(levels(as.factor(mapping[, 1, drop = TRUE])))
#     colors.light <- if(nGroups <= 8) {
#         suppressWarnings(RColorBrewer::brewer.pal(nGroups, name = "Set2"))[1:nGroups]
#         } else {
#             viridis::viridis(n = nGroups, alpha = 0.2)
#         }
#     colors.dark <- if(nGroups <= 8) {
#         suppressWarnings(RColorBrewer::brewer.pal(nGroups, name = "Dark2"))[1:nGroups]
#     } else {
#         viridis::viridis(n = nGroups, alpha = 1)
#     }
# } else {
#     colors.light <- viridis::viridis(n = nrow(mapping), alpha = 0.2)
#     colors.dark <- viridis::viridis(n = nrow(mapping), alpha = 1)
# }
# cat(paste("Your independent variable", colnames(mapping)[1], "has been detected as",
#           if(discrete) {"discrete."} else {"continuous."}
# ))


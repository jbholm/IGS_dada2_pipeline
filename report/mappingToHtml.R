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
      cat('\nTraceback:', file = stderr())
      cat('\n', file = stderr())
      traceback(2) # Print full traceback of function calls with all parameters. The 2 passed to traceback omits the outermost two function calls.
      if (!interactive()) quit(status = 1)
    },
    stringsAsFactors = FALSE
    )
invisible({
    initial.options <- commandArgs(trailingOnly = FALSE)
    wd <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
    fileConn <- file(file.path(wd, "config.txt"))
    lines <- readLines(fileConn)
    sapply(lines, function(line) {
        fields = strsplit(line, ": ")[[1]]
        if(fields[1] == "rlibpath") {
            .libPaths(c(fields[2], .libPaths()))
        }
    })

    lapply(initial.options, FUN = function(x) {
        cat(x)
        cat("\n")
    })

    suppressMessages(
        suppressWarnings(
            {
            require("optparse")
            require("kableExtra")
            require("stringi")
            library("magrittr")
            }
        )
    )
})
option_list = list(make_option(c("--map", "-m"), type="character", metavar = "MAPPING_FILE_PATH", help = "The mapping file, tab-delimited"),
                   make_option(c("--group_by", "-g"), type = "numeric", metavar = "COLUMN_INDEX", help = "Which column index to group by (optional)")
                   )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$map)) {
    stop("Mapping file must be provided using -m or --map")
}

mapping <- read.delim(opt$map, row.names=NULL, header=T, comment.char="", na.strings="")

colnames(mapping) = c("ID", "Barcode 1", "Barcode 2", "Description")



fileConn<-file("mapping.html.frag")
output <- mapping %>% kableExtra::kable() %>% kableExtra::kable_styling() %>% kableExtra::scroll_box(width = "100%", height = "100%")
output <- stringi::stri_replace_all(output, "<div class=\"kable\" style", regex = "^<div style")
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


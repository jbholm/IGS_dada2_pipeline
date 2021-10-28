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

config_file <- file.path(pipelineDir, "config.json")
config <- jsonlite::read_json(
    path = file.path(config_file)
)
.libPaths(config[["r-lib"]])

suppress_if_not_verbose <- function(block) {
    if (args$verbose) {
        block
    } else {
        suppressMessages(
            suppressWarnings(
                block
            )
        )
    }
}

project_meta <- function(runs = list(), samples = list()) {
    filepath <- ".meta.json"
    if (!file.exists(filepath)) {
        proj_info <- list()
    } else {
        proj_info <- jsonlite::read_json(
            path = filepath
        )
    }

    if (length(runs) > 0) {
        proj_info$runs[names(runs)] <- runs
    }
    if (length(samples) > 0) {
        proj_info$samples <- samples
    }
    if(length(runs) > 0 || length(samples) > 0) {
        jsonlite::write_json(proj_info, filepath, auto_unbox = T)
    }
    return(proj_info)
}


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

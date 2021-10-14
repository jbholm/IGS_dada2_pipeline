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
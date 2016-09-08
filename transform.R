## utility functions for changing values within a list or vector
## modified from transform.data.frame
transform.list <- function(`_data`,...) {
    e <- eval(substitute(list(...)), `_data`, parent.frame())
    tags <- names(e)
    inx <- match(tags, names(`_data`))
    matched <- !is.na(inx)
    if (any(matched)) {
        `_data`[inx[matched]] <- e[matched]
    }
    if (!all(matched)) 
        c(list(`_data`), e[!matched])
    else `_data`
}

transform.numeric <- function(x, ...) {
    transform(as.list(x))
}

stir_ER <- function(matrix, x, p, retcon = FALSE) {
    w <- x-1
    n <- ncol(matrix)

    ## Add zero because no self loops; diag(matrix) = 0
    before <- 1 * (runif(w) <= p)
    after <- 1 * (runif(n-x) <= p)
    matrix[x, ] = c(before, 0, after)

    if (retcon == TRUE) {
        stop("Retcon functionality missing.")
    }

    return(matrix)
}
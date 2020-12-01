grow_ER <- function(matrix, x, p, retcon = FALSE, directed = TRUE) {
    w <- x-1

    ## Add zero because no self loops; diag(matrix) = 0
    matrix[x, 1:x] = c(1 * (runif(w) <= p), 0)

    if (retcon == TRUE) {
        matrix[1:x, x] = c(1 * (runif(w) <= p), 0)
    } else {
        matrix[1:x, x] = 0
    }

    if (directed == FALSE) {
        matrix[1:x, x] = matrix[x, 1:x]
    }

    return(matrix)
}
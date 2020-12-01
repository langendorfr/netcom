stir_NM <- function(matrix, x, niches, connectance = 0.5, directed = TRUE) {
    ids <- (1:ncol(matrix))[-x]
    # w <- x-1

    # beta <- (1/connectance) - 1
    beta <- (1/(2*connectance)) - 1

    n_i <- niches[x]
    r_i <- 1-((1-runif(1))^(1/beta))

    r_i = r_i * n_i

    # if (r_i/2 > n_i) {
    #     r_i = 2 * n_i
    # }

    c_i <- runif(n = 1, min = r_i/2, max = n_i)

    range_min <- c_i - (r_i/2)
    range_max <- c_i + (r_i/2)

    interactions <- which( (niches[ids] >= range_min) & (niches[ids] <= range_max) )

    matrix[x, 1:x] = 0
    matrix[x, interactions] = 1

    matrix[1:x, x] = 0
    if (directed == FALSE) {
        matrix[1:x, x] = matrix[x, 1:x]
    }
  
    return(matrix)
}
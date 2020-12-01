grow_PA <- function(matrix, x, power, sum_v_max = "max", nascent_help = TRUE, retcon = FALSE, directed = TRUE) {
    w <- x-1

    out_degree <- rowSums(matrix[1:w, 1:w])
    in_degree <- colSums(matrix[1:w, 1:w])

    ## Give all nodes an extra edge so there is some nonzero probability of new nodes accumulating edges
    ## Note: This is likely less biased at small degrees than adding an edge only to zero degree nodes
    if (nascent_help == TRUE) {
        out_degree = out_degree + 1
        in_degree = in_degree + 1
    }

    out_power <- out_degree ^ power
    in_power <- in_degree ^ power

    if (sum_v_max == "sum") {
        out_ratio <- out_power / sum(out_power)
        in_ratio <- in_power / sum(in_power)
    } else if (sum_v_max == "max") {
        out_ratio <- out_power / max(out_power)
        in_ratio <- in_power / max(in_power)        
    } else {
        out_ratio <- 0
        in_ratio <- 0
    }


    out_new <- 1 * (runif(w) <= out_ratio)
    in_new <- 1 * (runif(w) <= in_ratio)

    matrix[x, 1:x] = c(in_new, 0)

    ## Add zero because no self loops; diag(matrix) = 0
    if (retcon == TRUE) {
        matrix[1:x, x] = c(out_new, 0)
    } else {
        matrix[1:x, x] = 0
    }

    if (directed == FALSE) {
        matrix[1:x, x] = matrix[x, 1:x]
    }

    return(matrix)
}
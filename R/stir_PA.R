stir_PA <- function(matrix, x, power, retcon = FALSE, sum_v_max = "max", nascent_help = TRUE) {
    # w <- x-1
    n <- ncol(matrix)

    out_degree <- rowSums(matrix[-x, -x])
    in_degree <- colSums(matrix[-x, -x])

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

    ## Prevent NAs when rows or columns are all zero
    if (any(is.na(out_ratio))) {
        out_ratio[which(is.na(out_ratio))] = 0
    }

    if (any(is.na(in_ratio))) {
        in_ratio[which(is.na(in_ratio))] = 0
    }

    out_new <- 1 * (runif(n-1) <= out_ratio)
    in_new <- 1 * (runif(n-1) <= in_ratio)

    if (retcon == TRUE) {
        stop("Retcon functionality missing.")
    }

    matrix[x, -x] = in_new

    return(matrix)
}
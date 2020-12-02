#' @title Stirs a Preferential Attachment Network
#'
#' @description Stirs an already existing network by rewiring a node according to the Preferential Attachment mechanism.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param x The ID of the node to be rewired (stirred).
#' 
#' @param power Power of attachment, which determines how much new nodes prefer to attach to nodes that have many edges compared to few. Needs to be positive.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#' 
#' @param sum_v_max Degree distributions must be normalized, either by their "max" or "sum". Defaults to "max".
#' 
#' @param nascent_help Should a single edge be added to the degree distribution of all nodes so that nodes with a zero in-degree can still have a chance of being attached to by new nodes. Defaults to TRUE.
#' 
#' @details Rewires a node in a network according to the Preferential Attachment mechanism.
#'
#' @return An adjacency matrix.
#' 
#' @references Barabási, A. L., & Albert, R. (1999). Emergence of scaling in random networks. science, 286(5439), 509-512.
#' 
#' @examples
#' size <- 10
#' existing_network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' new_network_prep <- matrix(0, nrow = size + 1, ncol = size + 1)
#' new_network_prep[1:size, 1:size] = existing_network
#' new_network <- stir_PA(matrix = new_network_prep, x = size + 1, power = 2.15)
#' 
#' @export

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
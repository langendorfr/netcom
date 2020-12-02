#' @title Grow a Preferential Attachment Network
#'
#' @description Grows an already existing network by adding a node according to the Preferential Attachment mechanism. Nodes can only attach to previously grown nodes.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param x The ID of the node to be grown.
#' 
#' @param power Power of attachment, which determines how much new nodes prefer to attach to nodes that have many edges compared to few. Needs to be positive.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Defaults to TRUE.
#'
#' @param sum_v_max Degree distributions must be normalized, either by their "max" or "sum". Defaults to "max".
#' 
#' @param nascent_help Should a single edge be added to the degree distribution of all nodes so that nodes with a zero in-degree can still have a chance of being attached to by new nodes. Defaults to TRUE.
#' 
#' @details Adds a node in a network according to the Preferential Attachment mechanism.
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
#' new_network <- grow_PA(matrix = new_network_prep, x = size + 1, power = 2.15)
#' 
#' @export

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
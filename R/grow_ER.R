#' @title Grow an Erdos-Renyi Random Network
#'
#' @description Grows an already existing network by adding a node according to the Erdos-Renyi random mechanism. Nodes can only attach to previously grown nodes.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param x The ID of the node to be grown.
#' 
#' @param p Probability possible edges exist. Needs to be between zero and one.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Defaults to TRUE.
#'
#' @details Different from Duplication & Mutation models in that edges can only be lost.
#'
#' @return An adjacency matrix.
#' 
#' @references Erdos, P. and Renyi, A., On random graphs, Publicationes Mathematicae 6, 290–297 (1959).
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' size <- 10
#' existing_network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' new_network_prep <- matrix(0, nrow = size + 1, ncol = size + 1)
#' new_network_prep[1:size, 1:size] = existing_network
#' new_network <- grow_ER(matrix = new_network_prep, x = size + 1, p = 0.5)
#' 
#' @export

grow_ER <- function(matrix, x, p, retcon = FALSE, directed = TRUE) {
    w <- x-1

    ## Add zero because no self loops; diag(matrix) = 0
    matrix[x, 1:x] = c(1 * (stats::runif(w) <= p), 0)

    if (retcon == TRUE) {
        matrix[1:x, x] = c(1 * (stats::runif(w) <= p), 0)
    } else {
        matrix[1:x, x] = 0
    }

    if (directed == FALSE) {
        matrix[1:x, x] = matrix[x, 1:x]
    }

    return(matrix)
}
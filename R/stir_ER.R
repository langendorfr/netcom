#' @title Stir an Erdos-Renyi Random Network
#'
#' @description Stirs an already existing network by rewiring a node according to the Erdos-Renyi random mechanism.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param x The ID of the node to be rewired (stirred).
#' 
#' @param p Probability possible edges exist. Needs to be between zero and one.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#'
#' @details Different from Duplication & Mutation models in that edges can only be lost.
#'
#' @return An adjacency matrix.
#' 
#' @references Erdos, P. and Renyi, A., On random graphs, Publicationes Mathematicae 6, 290–297 (1959).
#' 
#' @examples
#' size <- 10
#' existing_network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' new_network_prep <- matrix(0, nrow = size + 1, ncol = size + 1)
#' new_network_prep[1:size, 1:size] = existing_network
#' new_network <- stir_ER(matrix = new_network_prep, x = size + 1, p = 0.5)
#' 
#' @export

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
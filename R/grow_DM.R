#' @title Grow a Duplication and Mutation Network
#'
#' @description Grows an already existing network by adding a node according to the Duplication and Mutation mechanism. Nodes can only attach to previously grown nodes.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param x The ID of the node to be grown.
#' 
#' @param divergence Probability that the new node loses edges associated with the node it duplicates. Needs to be between zero and one.
#' 
#' @param mutation Probability that the new node gains edges not associated with the node it duplicates. Needs to be between zero and one.
#' 
#' @param link Probability that the new node attaches to the node it duplicates. Defaults to 0.
#' 
#' @param connected Binary argument determining if the newly grown node has to be connected to the existing network. Defaults to FALSE, to prevent rare computational slow-downs when it is unlikely to create a connected network. Defaults to FALSE.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Defaults to TRUE.
#'
#' @details Different from Duplication & Mutation models in that edges can only be lost.
#'
#' @return An adjacency matrix.
#' 
#' @references Ispolatov, I., Krapivsky, P. L., & Yuryev, A. (2005). Duplication-divergence model of protein interaction network. Physical review E, 71(6), 061911.
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' size <- 10
#' existing_network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' new_network_prep <- matrix(0, nrow = size + 1, ncol = size + 1)
#' new_network_prep[1:size, 1:size] = existing_network
#' new_network <- grow_DM(matrix = new_network_prep, x = size + 1, divergence = 0.5)
#' 
#' @export

grow_DM <- function(matrix, x, divergence, mutation = 0, link = 0, connected = FALSE, retcon = FALSE, directed = TRUE) {
    w <- x - 1
    DD <- function() {
        duplication <- sample(1:w, 1)
        matrix[, x] <- matrix[, duplication]
        matrix[x, ] <- matrix[duplication, ]

        ## Connect new node to copied node
        if (stats::runif(1) < link) {
        matrix[x, duplication] <- 1
        matrix[duplication, x] <- 1
        }

        for (i in 1:w) {

            if (matrix[x, i] == 1) {
                if (stats::runif(1) <= divergence) {
                    matrix[x, i] <- 0
                }
            } else if (matrix[x, i] == 0) {
                if (stats::runif(1) <= mutation) {
                    matrix[x, i] <- 0
                }
            } else {
                stop("Weighted edge detected. Only binary networks are supported in this release.")
            }

        }

        if (retcon == TRUE) {
            stop("Retcon functionality missing.")
        }

        if (directed == FALSE) {
            matrix[1:x, x] = matrix[x, 1:x]
        }

        return(matrix) ## DD
    }

    matrix = DD()

    ## Rerun until no entirely disconnected nodes
    if (connected == TRUE) {
        while (sum(matrix[, x]) == 0 && sum(matrix[x,]) == 0) {
            matrix = DD()
        }
    }

    return(matrix) ## grow_DD
}
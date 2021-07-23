#' @title Sitrs a Duplication and Divergence Network
#'
#' @description Stirs an already existing network by rewiring a node according to the Duplication and Divergence mechanism.
#'
#' @param matrix Existing network to be rewired (stirred).
#' 
#' @param x The ID of the node to be grown.
#' 
#' @param divergence Probability that the new node loses edges associated with the node it duplicates. Needs to be between zero and one.
#' 
#' @param link Probability that the new node attaches to the node it duplicates. Defaults to 0.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix.
#' 
#' @param force_connected Binary argument determining if the newly grown node has to be connected to the existing network. Defaults to FALSE, to prevent rare computational slow-downs when it is unlikely to create a connected network. Defaults to FALSE.
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
#' new_network <- stir_DD(matrix = new_network_prep, x = size + 1, divergence = 0.5)
#' 
#' @export

stir_DD <- function(matrix, x, divergence, link, directed, force_connected = FALSE) {
    ids <- (1:ncol(matrix))[-x]

    DD <- function() {
        duplication <- sample(ids, 1)
        matrix[, x] <- matrix[, duplication]
        matrix[x, ] <- matrix[duplication, ]

        ## Connect new node to copied node
        if (stats::runif(1) < link) {
        matrix[x,duplication] <- 1
        matrix[duplication,x] <- 1
        }

        for (i in ids) {
            ##!! Do this differently than grow_DD, which uses `if (matrix[i, x] == 1)`
            if (matrix[x, i] == 1) {
                if (stats::runif(1) < divergence) {
                    matrix[x, i] <- 0

                    if (!directed) {
                        matrix[i, x] <- 0
                    }   
                }
            }
        }
        
        return(matrix) ## DD
    }

    matrix = DD()

    if (force_connected) {
        ## Rerun until no entirely disconnected nodes
        ##!! Do this differently than grow_DD, which uses `while (sum(matrix[, x]) == 0 && sum(matrix[x,]) == 0) {`
        while (sum(matrix[x,]) == 0) {
        matrix = DD()
        }
    }

    return(matrix) ## grow_DD
}
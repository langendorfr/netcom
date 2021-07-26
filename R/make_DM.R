#' @title Make a Duplication and Mutation Network
#'
#' @description Make an already existing network according to the Duplication and Mutation mechanism.
#'
#' @param size Number of nodes in the network.
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list").
#' 
#' @param divergence Probability that the new node loses edges associated with the node it duplicates. Needs to be between zero and one.
#' 
#' @param mutation Probability that the new node gains edges not associated with the node it duplicates. Needs to be between zero and one.
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
#' # Network size (number of nodes)
#' size <- 10
#' 
#' # Divergence parameter
#' divergence <- 0.237
#' 
#' # Mutation parameter
#' mutation <- 0.1
#' 
#' # Make network according to the Duplication & Mutation mechanism
#' make_DM(size = size, net_kind = "matrix", divergence = divergence, mutation = mutation)
#' 
#' @export

make_DM <- function(size, net_kind, divergence, mutation, directed = FALSE) {
    if (net_kind == "matrix") {
        ## Start with pair of connected nodes
        matrix <- matrix(0, nrow = size, ncol = size)
        matrix[1,2] = 1
        matrix[2,1] = 1

        ## Start with node three because first two nodes are part of the initial network
        for (node in 3:size) {
            duplication <- sample(1:(node-1), 1)
            matrix[node, ] = matrix[duplication, ]

            if (directed == FALSE) {
                matrix[, node] = matrix[node, ]
            }

            ## Change each edge with probability divergence
            # edges <- which(matrix[node, ] != 0)
            edges <- matrix[node, 1:node]
            for (e in seq_along(edges)) {
                if (edges[e] == 1) {
                    if (stats::runif(1) <= divergence) {
                        matrix[node, e] = 0
                        
                        if (directed == FALSE) {
                            matrix[e, node] = 0
                        }
                    }
                } else if (edges[e] == 0) {
                    if (stats::runif(1) <= mutation) {
                        matrix[node, e] = 1

                        if (directed == FALSE) {
                            matrix[e, node] = 1
                        }
                    }
                } else {
                    stop("Weighted edge detected. Only binary networks are supported in this release.")
                }

            }
        }

        return(matrix)

    } else if (net_kind == "list") {

        stop("Edge lists are not supported for Duplication & Mutation models (make_DM) in this release. Set net_kind = `matrix`.")

    } else {
        stop("Unknown net_kind. Must be `list` or `matrix`.")
    }
}

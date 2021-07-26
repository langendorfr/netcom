#' @title Makes a Duplication and Divergence Network
#'
#' @description Makes a network according to the Duplication and Divergence mechanism.
#'
#' @param size Number of nodes in the network.
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list").
#' 
#' @param divergence Probability that the new node loses edges associated with the node it duplicates. Needs to be between zero and one.
#' 
#' @param directed Whether the target network is directed. Defaults to TRUE.
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
#' # Make network according to the Duplication & Divergence mechanism
#' make_DD(size = size, net_kind = "matrix", divergence = divergence)
#' 
#' @export

make_DD <- function(size, net_kind, divergence, directed = TRUE) {
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
            edges <- which(matrix[node, ] != 0)
            for (e in edges) {
                if (stats::runif(1) <= divergence) {
                    matrix[node, e] = 0

                    if (directed == FALSE) {
                        matrix[e, node] = 0
                    }
                }
            }
        }

        return(matrix)

    } else if (net_kind == "list") {
        edgelist <- matrix(nrow = 0, ncol = 2)
        edgelist = rbind(edgelist, c(1,2))
        edgelist = rbind(edgelist, c(2,1))

        ## Start with node three because first two nodes are part of the initial network
        for (node in 3:size) {
            duplication <- sample(1:(node-1), 1)
            edgelist_duplication <- matrix(edgelist[edgelist[,1] == duplication, ], ncol = 2)
            edgelist_duplication[,1] = node

            ## Change each edge with probability divergence
            for (e in 1:nrow(edgelist_duplication)) {
                if (stats::runif(1) <= divergence) {
                    edgelist_duplication = matrix(edgelist_duplication[-e,], ncol = 2)
                }
            }

            if (directed == FALSE) {
                edgelist_duplication_mirror <- cbind(edgelist_duplication[,2], edgelist_duplication[,1])
                edgelist_duplication = rbind(edgelist_duplication, edgelist_duplication_mirror)
            }

            edgelist = rbind(edgelist, edgelist_duplication)
        }

        return(edgelist)

    } else {
        stop("Unknown net_kind. Must be `list` or `matrix`.")
    }
}

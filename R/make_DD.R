#' @title 
#'
#' @description 
#'
#' @param 
#'
#' @details Different from Duplication & Mutation models in that edges can only be lost
#'
#' @return 
#' 
#' @references 
#' 
#' @examples
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
                if (runif(1) <= divergence) {
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
                if (runif(1) <= divergence) {
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

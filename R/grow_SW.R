#' @title Grow a Small-World Network
#'
#' @description Grows an already existing network by adding a node according to the Small-World mechanism. Nodes can only attach to previously grown nodes.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param x The ID of the node to be grown.
#' 
#' @param rewire Small-World parameter specifying the probability each edge is randomly rewired, allowing for the possiblity of bridges between connected communities.
#' 
#' @param connected Binary argument determining if the newly grown node has to be connected to the existing network. Defaults to FALSE, to prevent rare computational slow-downs when it is unlikely to create a connected network. Defaults to False.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Defaults to TRUE.
#' 
#' @details Grows a node in a network according to the Small-World mechanism.
#'
#' @return An adjacency matrix.
#' 
#' @references Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of ‘small-world’networks. nature, 393(6684), 440-442.
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' size <- 10
#' existing_network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' new_network_prep <- matrix(0, nrow = size + 1, ncol = size + 1)
#' new_network_prep[1:size, 1:size] = existing_network
#' new_network <- grow_SW(matrix = new_network_prep, x = size + 1, rewire = 0.213)
#' 
#' @export

grow_SW <- function(matrix, x, rewire, connected = FALSE, retcon = FALSE, directed = TRUE) {
    w <- x-1

    out_degree <- rowSums(as.matrix(matrix[1:w, 1:w]))
    in_degree <- colSums(as.matrix(matrix[1:w, 1:w]))

    out_avg <- mean(out_degree) %>% round(digits = 0)
    in_avg <- mean(in_degree) %>% round(digits = 0)

    out_new <- 1 * ((out_avg - out_degree) >= 1)
    in_new <- 1 * ((in_avg - in_degree) >= 1)

    ## Do not allow no links
    if (connected == TRUE) {
        if (max(out_new) == 0) {
            out_id <- sample(seq_along(out_new), size = 1)
            in_id <- sample(seq_along(in_new), size = 1)

            out_new[out_id] = 1
            in_new[in_id] = 1
        }
    }

    ## Add zero because no self loops; diag(matrix) = 0
    matrix[x, 1:x] = c(out_new, 0)
    
    if (retcon == TRUE) {
        matrix[1:x, x] = c(in_new, 0)
    } else {
        matrix[1:x, x] = 0
    }

    for (col in 1:w) {
    
        ## Do not attempt rewiring when all edges exist
        out_edges <- which(matrix[x, 1:w] != 0)
        if (length(out_edges) != w) {
            for (e in out_edges) {
                if (stats::runif(1) <= rewire) {
                    e_possible <- which(matrix[x, 1:w] == 0)
                    e_new <- sample(e_possible, 1)
                    
                    matrix[x, e] = 0
                    matrix[x, e_new] = 1
                }
            }
        }
    }

    if (directed == FALSE) {
        matrix[1:x, x] = matrix[x, 1:x]
    }

    return(matrix)
}
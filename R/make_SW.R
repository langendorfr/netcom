#' @title Makes a Small-World Network
#'
#' @description Make an already existing network according to the Small-World mechanism.
#'
#' @param size The number of nodes in the network. Must be a positive integer.
#' 
#' @param net_kind The format of the network. Currently must be either `matrix` or `list`.x
#' 
#' @param rewire Small-World parameter specifying the probability each edge is randomly rewired, allowing for the possiblity of bridges between connected communities.
#' 
#' @param neighborhood The range of nodes that form connected communities. Note: This implementation results in overlap of communities.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Defaults to TRUE.
#' 
#' @details Rewires a node in a network according to the Small-World mechanism.
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
#' new_network <- make_SW(matrix = new_network_prep, x = size + 1, rewire = 0.213)
#' 
#' @export

make_SW <- function(size, net_kind, rewire, neighborhood, directed = FALSE) {
    ## Default neighborhood is one-tenth the network's size
    if (missing(neighborhood)) {
        neighborhood = max(1, round(size/10))
    }
    
    ## Check parameter viability
    if ((rewire < 0) | (rewire > 1)) {
        stop("The rewire parameter must be between 0 and 1.")
    }

    if ((neighborhood < 1) | (neighborhood > size)) {
        stop("The neighborhood parameter must be between 1 and the `size` of the network.")
    }
    
    if (net_kind == "matrix") {
        matrix <- matrix(0, nrow = size, ncol = size)

        for (node in 1:size) {
            ids <- seq(from = node - neighborhood, to = node + neighborhood, by = 1)
            
            ## Wrap ids so node = 1 and node = size are adjacent
            ids[ids < 1] = size - abs(ids[ids < 1])
            ids[ids > size] = ids[ids > size] - size

            matrix[node, ids] = 1
        }

        ## Remove self-loops
        diag(matrix) = 0

        ## Rewiring, which is where `directed` comes into play
        edge_ids <- which(matrix != 0, arr.ind = TRUE)
        
        if (directed == TRUE) {
            for (edge in 1:nrow(edge_ids)) {
                if (stats::runif(1) <= rewire) {
                    ## Prevent rewiring to self or current target
                    possible_ids <- which(matrix[edge_ids[edge, 1], ] == 0)
                    impossible_ids <- edge_ids[edge,]
                    possible_ids = possible_ids[-which(possible_ids %in% impossible_ids)]

                    ## If all edges already exist no rewiring is possible
                    if (length(possible_ids) > 0) {
                        new_target <- sample(x = possible_ids, size = 1)

                        matrix[edge_ids[edge, 1], edge_ids[edge, 2]] = 0
                        matrix[edge_ids[edge, 1], new_target] = 1
                    }
                }
            }

        } else if (directed == FALSE) {
            for (edge in 1:nrow(edge_ids)) {
                ## The starting regular network is symmetric across the diagonal
                ## Only consider the lower diagonal and apply changes to both edges
                ## The is asymmetric in which node is kept in the edge, but resulting network structures are invariant to this choice
                if (edge_ids[edge, 2] < edge_ids[edge, 1]) {
                    if (stats::runif(1) <= rewire) {
                        ## Prevent rewiring to self or current target
                        possible_ids <- which(matrix[edge_ids[edge, 1], ] == 0)
                        impossible_ids <- edge_ids[edge,]
                        possible_ids = possible_ids[-which(possible_ids %in% impossible_ids)]
                        
                        ## If all edges already exist no rewiring is possible
                        if (length(possible_ids) > 0) {
                            new_target <- sample(x = possible_ids, size = 1)

                            matrix[edge_ids[edge, 1], edge_ids[edge, 2]] = 0
                            matrix[edge_ids[edge, 2], edge_ids[edge, 1]] = 0

                            matrix[edge_ids[edge, 1], new_target] = 1
                            matrix[new_target, edge_ids[edge, 1]] = 1
                        }
                    }
                }
            }

        } else {
            stop("The parameter `directed` is a logical value that must be either TRUE or FALSE.")
        }

        return(matrix)

    } else if (net_kind == "list") {
        stop("Currently only net_kind = `matrix` is supported.")

    } else {
        stop("Unknown net_kind. Must be `list` or `matrix`.")
    }
}

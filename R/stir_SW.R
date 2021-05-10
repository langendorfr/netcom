#' @title Stirs a Small-World Network
#'
#' @description Stirs an already existing network by rewiring a node according to the Small-World mechanism.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param x The ID of the node to be grown.
#' 
#' @param rewire Small-World parameter specifying the probability each edge is randomly rewired, allowing for the possiblity of bridges between connected communities.
#' 
#' @details Rewires a node in a network according to the Small-World mechanism.
#'
#' @return An adjacency matrix.
#' 
#' @references Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of ‘small-world’networks. nature, 393(6684), 440-442.
#' 
#' @examples
#' size <- 10
#' existing_network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' new_network_prep <- matrix(0, nrow = size + 1, ncol = size + 1)
#' new_network_prep[1:size, 1:size] = existing_network
#' new_network <- stir_SW(matrix = new_network_prep, x = size + 1, rewire = 0.213)
#' 
#' @export

stir_SW <- function(matrix, x, rewire) {
    ids <- (1:ncol(matrix))[-x]

    out_degree <- rowSums(as.matrix(matrix[-x, -x]))
    # in_degree <- colSums(as.matrix(matrix[-x, -x]))

    out_avg <- mean(out_degree) %>% round(digits = 0)
    # in_avg <- mean(in_degree) %>% round(digits = 0)

    out_new <- 1 * ((out_avg - out_degree) >= 1)
    # in_new <- 1 * ((in_avg - in_degree) >= 1)

    ## Do not allow no links
    if (max(out_new) == 0) {
        out_id <- sample(seq_along(out_new), size = 1)
        # in_id <- sample(seq_along(in_new), size = 1)

        out_new[out_id] = 1
        # in_new[in_id] = 1        
    }

    ## Add zero because no self loops; diag(matrix) = 0
    matrix[x, -x] = out_new
    # matrix[x, 1:x] = c(in_new, 0)

    # out_degree_grown <- rowSums(matrix[1:x, 1:x])

    # for (row in 1:x) {
    out_edges <- which(matrix[x, ] != 0)
    if (x %in% out_edges) {
        self_id <- which(out_edges == x)
        out_edges = out_edges[-self_id]
    }

    ## -1 because no self loops so possible number of out_edges is ncol(matrix)-1
    if ( length(out_edges) != (ncol(matrix)-1) ) {
        for (e in out_edges) {
            
            ## Do not attempt rewiring when all edges exist
            ## Do not allow self-loops so remove column x
            e_possible <- which(matrix[x, -x] == 0)

            if (runif(1) <= rewire && length(e_possible) != 0) {
                e_new <- sample(e_possible, 1)
                matrix[x, e] = 0
                matrix[x, e_new] = 1
            }
        }
    }


    return(matrix)
}
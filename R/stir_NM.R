#' @title Stirs a Niche Model Network
#'
#' @description Stirs an already existing network by rewiring a node according to the Niche Model mechanism.
#'
#' @param matrix Existing network to experience rewiring (stirring).
#' 
#' @param x The ID of the node to be grown.
#' 
#' @param niches Vector of length x, with values between zero and one corresponding to each node's niche.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Defaults to TRUE.
#' 
#' @param connectance Niche Model parameter specifying the expected connectivity of the network, which determines for a given node the niche space window within which it attaches to every other node. Defaults to 0.2.
#'
#' @details Stirs a node in a Niche Model network.
#'
#' @return An adjacency matrix.
#' 
#' @references Williams, R. J., & Martinez, N. D. (2000). Simple rules yield complex food webs. Nature, 404(6774), 180-183.
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' size <- 10
#' existing_network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' new_network_prep <- matrix(0, nrow = size + 1, ncol = size + 1)
#' new_network_prep[1:size, 1:size] = existing_network
#' new_network <- stir_NM(matrix = new_network_prep, x = size + 1, power = 2.15)
#' 
#' @export

stir_NM <- function(matrix, x, niches, directed, connectance = 0.2) {
    ids <- (1:ncol(matrix))[-x]
    # w <- x-1

    # beta <- (1/connectance) - 1
    beta <- (1/(2*connectance)) - 1

    n_i <- niches[x]
    r_i <- 1-((1-stats::runif(1))^(1/beta))

    r_i = r_i * n_i

    # if (r_i/2 > n_i) {
    #     r_i = 2 * n_i
    # }

    c_i <- stats::runif(n = 1, min = r_i/2, max = n_i)

    range_min <- c_i - (r_i/2)
    range_max <- c_i + (r_i/2)

    interactions <- which( (niches[ids] >= range_min) & (niches[ids] <= range_max) )

    matrix[x, 1:x] = 0
    matrix[x, interactions] = 1

    matrix[1:x, x] = 0
    if (!directed) {
        matrix[1:x, x] = matrix[x, 1:x]
    }
  
    return(matrix)
}

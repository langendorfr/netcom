#' @title Make a Niche Model network
#'
#' @description Creates a single network according to the Niche Model. Can be directed or undirected, but is always unweighted.
#'
#' @param size The number of nodes in the network. Must be a positive integer.
#' 
#' @param niches A vector of numbers specifying the niche of each member of the system (node). Each niche value must be element of [0,1].
#' 
#' @param connectance Defaults to 0.5. The ratio of actual interactions to possible interactions. Effects the beta distributed width of niche values each member of the system (node) interacts with.
#'
#' @param directed Defaults to TRUE. If FALSE all interactions will be made symmetric. Note that the process of creating interactions is unaffected by this choice.
#' 
#' @return An interaction matrix format of a Niche Model network.
#' 
#' @references Williams, R. J., & Martinez, N. D. (2000). Simple rules yield complex food webs. Nature, 404(6774), 180-183.
#'
#' @examples
#' # Network size (number of nodes)
#' size <- 10
#' 
#' # Create niche values for each member of the system (node)
#' niches <- runif(n = size)
#' 
#' # Make network according to the Niche Model
#' make_NM(size = size, niches = niches)
#' 
#' @export

make_NM <- function(size, niches, connectance = 0.5, directed = TRUE) {
    matrix <- matrix(0, 
                     nrow = size, 
                     ncol = size)

    beta <- (1/connectance) - 1

    for (x in 1:size) {
        n_i <- niches[x]
        r_i <- 1-((1-runif(1))^(1/beta))

        if (r_i/2 > n_i) {
            r_i = 2 * n_i
        }

        c_i <- runif(n = 1, min = r_i/2, max = n_i)
        
        range_min <- c_i - (r_i/2)
        range_max <- c_i + (r_i/2)

        interactions <- which( (niches >= range_min) & (niches <= range_max) )

        matrix[x,] = 0
        matrix[x, interactions] = 1

        if (directed == FALSE) {
            matrix[, x] = matrix[x,]
        }
    } ## for (x in 1:size)

    return(matrix)
}

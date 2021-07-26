#' @title Make a Niche Model network
#'
#' @description Creates a single network according to the Niche Model. Can be directed or undirected, but is always unweighted.
#'
#' @param size The number of nodes in the network. Must be a positive integer.
#' 
#' @param net_kind The format of the network. Currently must be either `matrix` or `list`.
#' 
#' @param niches A vector of numbers specifying the niche of each member of the system (node). Each niche value must be element of [0,1].
#' 
#' @param connectance Defaults to 0.5. The ratio of actual interactions to possible interactions. Effects the beta distributed width of niche values each member of the system (node) interacts with.
#'
#' @param directed If FALSE all interactions will be made symmetric. Note that the process of creating interactions is unaffected by this choice. Defaults to TRUE.
#' 
#' @param grow Binary argument that determines if the network should be made in a growing fashion, where nodes' edges are added in order of their niches and can only attach to previously considered nodes. Defaults to FALSE.
#' 
#' @return An interaction matrix format of a Niche Model network.
#' 
#' @references Williams, R. J., & Martinez, N. D. (2000). Simple rules yield complex food webs. Nature, 404(6774), 180-183.
#'
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' # Network size (number of nodes)
#' size <- 10
#' 
#' # Create niche values for each member of the system (node)
#' niches <- stats::runif(n = size)
#' 
#' # Make network according to the Niche Model
#' make_NM(size = size, niches = niches)
#' 
#' @export

make_NM <- function(size, niches, net_kind = "matrix", connectance = 0.1, directed = TRUE, grow = FALSE) {

    if (net_kind == "matrix") {
        matrix <- matrix(0, 
                         nrow = size,
                         ncol = size)
    } else if (net_kind == "list") {
        edgelist <- matrix(nrow = 0,
                           ncol = 2)
    } else {
        stop("Unknown net_kind. Must be `list` or `matrix`.")
    }

    # beta <- (1/connectance) - 1
    beta <- (1/(2*connectance)) - 1

    for (x in 1:size) {
        n_i <- niches[x]
        
        r_i <- 1-((1-stats::runif(1))^(1/beta))

        ## NEW
        r_i = r_i * n_i



        # if (r_i/2 > n_i) {
        #     r_i = 2 * n_i
        # }

        c_i <- stats::runif(n = 1, 
                     min = min(n_i, r_i/2), ## Do not allow c_i to be greater than n_i
                     max = n_i)

        range_min <- c_i - (r_i/2)
        range_max <- c_i + (r_i/2)

        interactions <- which( (niches >= range_min) & (niches <= range_max) )

        if (grow) {
            interactions = interactions[which(interactions <= x)]
        }


        if (net_kind == "matrix") {
            matrix[x,] = 0
            matrix[x, interactions] = 1

            if (directed == FALSE) {
                matrix[, x] = matrix[x,]
            }

        } else if (net_kind == "list") {
            for (i in interactions) {
                edgelist = rbind(edgelist, c(x, i))

                if (directed == FALSE) {
                    edgelist = rbind(edgelist, c(i, x))
                }
            }

        } else {
            stop("Unknown net_kind. Must be `list` or `matrix`.")
        }

    } ## for (x in 1:size)

    ## Output
    if (net_kind == "matrix") {
        return(matrix)

    } else if (net_kind == "list") {
        return(edgelist)

    } else {
        stop("Unknown network kind. Must be `list` or `matrix`.")
    }

}

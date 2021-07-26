#' @title Make a Mixture Mechanism Network
#'
#' @description Creates a network by iteratively adding or rewiring nodes, each capable of attaching to existing nodes according to a user-specified mechanism.
#' 
#' @param mechanism A vector of mechanism names corresponding to the mechanisms each node acts in accordance with. Note that the first two mechanisms are irrelevant because the first two nodes default to connecting to each other. Currently supported mechanisms: "ER" (Erdos-Renyi random), "PA", (Preferential Attachment), "DD", (Duplication and Divergence), "DM" (Duplication and Mutation), "SW", (Small-World), and "NM" (Niche Model).
#' 
#' @param directed A binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Either a single value or a vector of values the same length as the mechanism input vector.
#' 
#' @param parameter Parameter of each node's mechanism. Either a single value or a vector of values the same length as the mechanism input vector.
#' 
#' @param kind Either `grow` or `rewire`, and determines if the nodes specified in the mechanism input vector are to be rewired or grown. Either a single value or a vector of values the same length as the mechanism input vector. The number of `grow` nodes, excluding the first two which are always a pair of bidirectionally connected nodes, is the size of the final network.
#'
#' @param size Typically not specified. The size of the network depends on how many `grow` events are part of the `kind` input sequence. This should only be used when all four components of the network evolution (`mechanism`, `kind`, `parameter`, and `directed`) are single name inputs instead of vectors.
#' 
#' @param niches Used by the Niche Model to determine which nodes interact. Needs to be a vector of the same length as the number of nodes, and range between zero and one.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#' 
#' @param link_DD Defaults to 0. A second parameter in the DD (Duplication & Divergence). Currently only one parameter per mechanism can be specified.
#' 
#' @param link_DM Defaults to 0. A second parameter in the DM (Duplication & Mutation). Currently only one parameter per mechanism can be specified.
#' 
#' @param force_connected Defaults to FALSE. Determines if nodes can be added to the growing network that are disconnected. If TRUE, this is prevented by re-determining the offending node's edges until the network is connected.
#' 
#' @details This function grows, one node at a time, a mixture mechanism network. As each node is added to the growing network it can attach to existing nodes by its own node-specific mechanism. A sequence of mechanism names must be provided. Note: Currently each mechanism is assumed to have a single governing parameter.
#'
#' @return An unweighted mixture mechanism adjacency matrix.
#' 
#' @references Langendorf, R. E., & Burgess, M. G. (2020). Empirically Classifying Network Mechanisms. arXiv preprint arXiv:2012.15863.
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' # Start by creating a sequence of network evolutions. 
#' # There are four components to this sequence that can each be defined for every step 
#' # in the network's evolution. Or, you can also specify a component once which will 
#' # be used for every step in the newtwork's evolution.
#' 
#' mechanism <- c(
#'     rep("ER", 7),
#'     rep("PA", 2),
#'     rep("ER", 3)
#' )
#' 
#' kind <- c(
#'     rep("grow", 7),
#'     rep("rewire", 2),
#'     rep("grow", 3)
#' )
#' 
#' parameter <- c(
#'     rep(0.3, 7),
#'     rep(2, 2),
#'     rep(0.3, 3)
#' )

#' directed <- c(
#'     rep(TRUE, 7),
#'     rep(FALSE, 2),
#'     rep(TRUE, 3)
#' )
#' 
#' # Simulate a network according to the rules of this system evolution.
#' network <- make_Mixture(
#'      mechanism = mechanism, 
#'      kind = kind, 
#'      parameter = parameter, 
#'      directed = directed
#' )
#' 
#' @export

make_Mixture <- function(mechanism, directed, parameter, kind, size, niches, retcon = FALSE, link_DD = 0, link_DM = 0, force_connected = FALSE) {

    ## Check that the input components are of matching lengths if not a single value
    component_sizes <- {}
    if (!missing(mechanism)) {
        component_sizes = c(component_sizes, length(mechanism))
    }
    if (!missing(directed)) {
        component_sizes = c(component_sizes, length(directed))
    }
    if (!missing(parameter)) {
        component_sizes = c(component_sizes, length(parameter))
    }
    if (!missing(kind)) {
        component_sizes = c(component_sizes, length(kind))
    }

    if (length(component_sizes) == 0) {
            stop("All four input components (`mechanism`, `kind`, `parameter`, and `directed`) are missing. Specify at least one of them.")
    }

    system_sizes <- unique(component_sizes)
    system_sizes_vectors <- system_sizes[which(system_sizes != 1)]

    if (length(system_sizes_vectors) > 1) {
        stop("Make sure all four input components are either a single value or vectors of the same length. These correspond to node events in the evolution of the network. Every event must have all four components.")
    }

    if (max(system_sizes) > 1) {
        ## Okay to overwrite the input variable `size` because the previous if statement
        ## ensures that either `kind` is a single value or the same length as any other
        ## input component.
        size = max(system_sizes)
    }

    if (length(mechanism) == 1) {
        if (missing(size)) {
            stop("If all four input components (`mechanism`, `kind`, `parameter`, and `directed`) are a single value, use the `size` input to specify the size of the network.")
        }
        mechanism = rep(mechanism, size)
    }

    ## Handle missing kind which can be assumed to be "grow"
    ## All other information (mechanism, directed, and kind) should be specified
    if (missing(kind)) {
        kind <- rep("grow", length(mechanism))
    }

    ## Handle constant inputs
    if (length(parameter) == 1) {
        parameter = rep(parameter, length(mechanism))
    }
    if (length(kind) == 1) {
        kind = rep(kind, length(mechanism))
    }
    if (length(directed) == 1) {
        directed = rep(directed, length(mechanism))
    }

    ## Okay to overwrite the input variable `size` (possible for the second
    ## time) because the network evolution has been specified or an error
    ## thrown by this point.
    # if (size != sum(kind == "grow")) {
    #     stop("The specified `size` is different from the specified number of `grow` events.")
    # }
    size = sum(kind == "grow")
    matrix <- matrix(0, size, size)

    ## Start with the first two nodes connected
    matrix[1,2] = 1
    matrix[2,1] = 1

    ## Assign niches for NM (the Niche Model) so they will be constant across calls to grow_CM
    if (missing(niches)) {
        niches <- stats::runif(size) %>% sort()
    } else {
        # Sort to be safe with use-supplied niches
        niches = niches %>% sort()
    }

    ## Keep track of growing size, which is different than `x` because rewiring events do not add nodes
    s <- 2

    for (x in 3:length(mechanism)) {
        if (kind[x] == "grow") {
            s = s + 1
        } else if (kind[x] == "rewire") {
            x_rewire <- sample(1:s, 1)
        } else {
            stop("ERROR: Kind of network evolution not specified correctly. `kind` must be a vector of either `grow` or `rewire`, for growing and rewiring events respectively.")
        }

        if (kind[x] == "grow") {
            if (mechanism[x] == "ER") {
                matrix = grow_ER(matrix, s, p = parameter[x], retcon = retcon, directed = directed[x]) #TRUE)
            } else if (mechanism[x] == "PA") {
                matrix = grow_PA(matrix, s, power = parameter[x], retcon = retcon, directed = directed[x]) #TRUE)
            } else if (mechanism[x] == "DD") {
                matrix = grow_DD(matrix, s, divergence = parameter[x], link = link_DD, directed = directed[x]) #FALSE)
            } else if (mechanism[x] == "DM") {
                matrix = grow_DM(matrix, s, divergence = parameter[x], mutation = parameter[x], link = link_DM, directed = directed[x]) #FALSE)
            } else if (mechanism[x] == "SW") {
                matrix = grow_SW(matrix, s, rewire = parameter[x], retcon = retcon, directed = directed[x]) #FALSE)
            } else if (mechanism[x] == "NM") {
                matrix = grow_NM(matrix, s, connectance = parameter[x], niches = niches, retcon = retcon, directed = directed[x]) #TRUE)
            } else {
                stop("ERROR: Model not specified.")
            }
        } else if (kind[x] == "rewire") {
            if (mechanism[x] == "ER") {
                matrix = stir_ER(matrix = matrix, x = x_rewire, p = parameter[x], directed = directed[x])
            } else if (mechanism[x] == "PA") {
                matrix = stir_PA(matrix = matrix, x = x_rewire, power = parameter[x], directed = directed[x])
            } else if (mechanism[x] == "DD") {
                matrix = stir_DD(matrix = matrix, x = x_rewire, divergence = parameter[x], link = link_DD, force_connected = force_connected, directed = directed[x])
            } else if (mechanism[x] == "DM") {
                matrix = stir_DM(matrix = matrix, x = x_rewire, divergence = parameter[x], mutation = parameter[x], link = link_DM, force_connected = force_connected, directed = directed[x])
            } else if (mechanism[x] == "NM") {
                matrix = stir_NM(matrix = matrix, x = x_rewire, niches = niches, connectance = parameter[x], directed = directed[x])
            } else if (mechanism[x] == "SW") {
                matrix = stir_SW(matrix = matrix, x = x_rewire, rewire = parameter[x], directed = directed[x]   )
            } else {
                stop("ERROR: Model not specified.")
            }
        }
    }

    return(matrix)
}

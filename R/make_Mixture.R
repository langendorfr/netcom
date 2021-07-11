#' @title Make a Mixture Mechanism Network
#'
#' @description Creates a network by iteratively adding or rewiring nodes, each capable of attaching to existing nodes according to a user-specified mechanism.
#'
#' @param sequence A vector of mechanism names corresponding to the mechanisms each node acts in accordance with. Needs to be the same length as the number of nodes in the network. Note that the first two mechanisms are irrelevant because the first two nodes default to connecting to each other. Currently supported mechanisms: "ER" (Erdos-Renyi random), "PA", (Preferential Attachment), "DD", (Duplication and Divergence), "DM" (Duplication and Mutation), "SW", (Small-World), and "NM" (Niche Model).
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix.
#'
#' @param niches Used by the Niche Model to determine which nodes interact. Needs to be a vector of the same length as the number of nodes, and range between zero and one.
#' 
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
#' 
#' @details This function grows, one node at a time, a mixture mechanism network. As each node is added to the growing network it can attach to existing nodes by its own node-specific mechanism. A sequence of mechanism names must be provided. Note: Currently each mechanism is assumed to have a single governing parameter.
#'
#' @return An unweighted mixture mechanism adjacency matrix.
#' 
#' @references Langendorf, R. E., & Burgess, M. G. (2020). Empirically Classifying Network Mechanisms. arXiv preprint arXiv:2012.15863.
#' 
#' @examples
#' # Mechanisms
#' sequence <- c("gER", "gER", "gER", "rPA", "rPA", "gER", "gER")
#' network <- make_Mixture(sequence)
#' 
#' @export

make_Mixture <- function(mechanism, directed, parameter, kind, niches, retcon = FALSE, link_DD = 0, link_DM = 0, force_connected = FALSE) {

    ## Old way of all information being in a single `sequence` vector that gets parsed
    # if (ncol(as.data.frame(do.call(rbind, strsplit(sequence, "_")))) == 3) {
    #     sequence_pars <- as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[2]] %>% as.numeric()
    #     if (all(as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[3]] %in% c("grow", "rewire"))) {
    #         sequence = paste0(substr(as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[3]], 1, 1), as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[1]])
    #     } else {
    #         stop("Sequence's third variable, after the second underscore, must be either 'grow' or 'rewire'. For example: sequence <- c('SW_0.15_grow', 'SW_0.6_rewire', 'SW_0.2_grow', 'SW_0.2_grow')")
    #     }
    # } else if (ncol(as.data.frame(do.call(rbind, strsplit(sequence, "_")))) == 2) {
    #     sequence_pars <- as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[2]] %>% as.numeric()
    #     sequence = paste0("g", as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[1]])
    # }

    ## New way of having a vector for each piece of information
    if (missing(directed)) {
        directed <- rep(TRUE, length(mechanism))
    }

    if (missing(kind)) {
        kind <- rep("grow", length(mechanism))
    }

    ## Handle constant parameter
    if (length(parameter) == 1) {
        parameters = rep(parameter, length(mechanism))
    }

    size <- sum(kind == "grow")
    matrix <- matrix(0, size, size)

    ## Start with the first two nodes connected
    matrix[1,2] = 1
    matrix[2,1] = 1

    ## Assign niches for NM (the Niche Model) so they will be constant across calls to grow_CM
    if (missing(niches)) {
        niches <- runif(size) %>% sort()
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

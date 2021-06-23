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

make_Mixture <- function(sequence, directed niches, mechanism_kind = "grow", retcon = FALSE, link_DD = 0, link_DM = 0) {

    # if (ncol(as.data.frame(do.call(rbind, strsplit(sequence, "_")))) != 1) {
    sequence_pars <- as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[2]] %>% as.numeric()
    sequence = paste0(as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[3]], as.data.frame(do.call(rbind, strsplit(sequence, "_")))[[1]])
    # }

    size <- sum(substr(sequence, 1, 1) == "g")
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

    for (x in 3:length(sequence)) {
        kind <- substr(sequence[x], 1, 1)
        if (kind == "g") {
            s = s + 1
        } else if (kind == "r") {
            x_rewire <- sample(1:s, 1)
        } else {
            stop("ERROR: Kind of network evolution not specified correctly. `sequence` must begin with either `g` or `r`, for growing and rewiring events respectively.")
        }

        if (sequence[x] == "gER") {
            matrix = grow_ER(matrix, s, p = sequence_pars[x], retcon = retcon, directed = directed) #TRUE)
        } else if (sequence[x] == "gPA") {
            # print(paste0("asldkfjasldkfjasdlfj:    ", sequence_pars[x]))
            matrix = grow_PA(matrix, s, power = sequence_pars[x], retcon = retcon, directed = directed) #TRUE)
        } else if (sequence[x] == "gDD") {
            matrix = grow_DD(matrix, s, divergence = sequence_pars[x], link = link_DD, directed = directed) #FALSE)
        } else if (sequence[x] == "gDM") {
            matrix = grow_DM(matrix, s, divergence = sequence_pars[x], mutation = sequence_pars[x], link = link_DM, directed = directed) #FALSE)
        } else if (sequence[x] == "gSW") {
            matrix = grow_SW(matrix, s, rewire = sequence_pars[x], retcon = retcon, directed = directed) #FALSE)
        } else if (sequence[x] == "gNM") {
            matrix = grow_NM(matrix, s, connectance = sequence_pars[x], niches = niches, retcon = retcon, directed = directed) #TRUE)
        } else if (sequence[x] == "rER") {
            matrix = stir_ER(matrix = matrix, x = x_rewire, p = sequence_pars[x], directed = directed)
        } else if (sequence[x] == "rPA") {
            matrix = stir_PA(matrix = matrix, x = x_rewire, power = sequence_pars[x], directed = directed)
        } else if (sequence[x] == "rDD") {
            matrix = stir_DD(matrix = matrix, x = x_rewire, divergence = sequence_pars[x], link = link_DD, force_connected = force_connected, directed = directed)
        } else if (sequence[x] == "rDM") {
            matrix = stir_DM(matrix = matrix, x = x_rewire, divergence = sequence_pars[x], mutation = sequence_pars[x], link = link_DM, force_connected = force_connected, directed = directed)
        } else if (sequence[x] == "rNM") {
            matrix = stir_NM(matrix = matrix, x = x_rewire, niches = niches, connectance = sequence_pars[x], directed = directed)
        } else if (sequence[x] == "rSW") {
            matrix = stir_SW(matrix = matrix, x = x_rewire, rewire = sequence_pars[x], directed = directed)
        } else {
            stop("ERROR: Model not specified.")
        }
    }

    return(matrix)
}

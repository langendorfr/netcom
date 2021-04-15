#' @title Grow a Mixture Mechanism Network
#'
#' @description Creates a network by iteratively adding nodes, each capable of attaching to existing nodes according to a user-specified mechanism.
#'
#' @param matrix Existing network to experience growth.
#' 
#' @param sequence A vector of mechanism names corresponding to the mechanisms each node acts in accordance with. Needs to be the same length as the number of nodes in the network. Note that the first two mechanisms are irrelevant because the first two nodes default to connecting to each other. Currently supported mechanisms: "ER" (Erdos-Renyi random), "PA", (Preferential Attachment), "DD", (Duplication and Divergence), "DM" (Duplication and Mutation), "SW", (Small-World), and "NM" (Niche Model).
#' 
#' @param stirs Number of times to stir every node in the network. This should be more than one because the order matters. Stirring many times makes the order matter less.
#' 
#' @param p_ER Erdos-Renyi parameter specifying the probability each possible edge actually exists. Defaults to 0.5.
#' 
#' @param power_PA Preferential Attachment parameter specifying the power of attachment, which determines how much new nodes prefer to attach to nodes that are already attached to by other nodes. Defaults to 2.
#' 
#' @param divergence_DD Duplication and Divergence parameter specifying the probability of new nodes losing edges that exist in the node they duplicated. Defaults to 0.1.
#' 
#' @param link_DD Duplication and Divergence parameter specifying the probability of an edge between a new node and the node it duplicated.
#' 
#' @param divergence_DM Duplication and Mutation parameter specifying the probability of new nodes losing edges that exist in the node they duplicated. Defaults to 0.
#' 
#' @param mutation_DM Duplication and Mutation parameter specifying the probability of new nodes gaining edges that did not exist in the node they duplicated. Defaults to 0.1.
#' 
#' @param link_DM Duplication and Mutation parameter specifying the probability of an edge between a new node and the node it duplicated. Defaults to 0.
#' 
#' @param rewire_SW Small-World parameter specifying the probability each edge is randomly rewired, allowing for the possiblity of bridges between connected communities. Defaults to 0.1.
#' 
#' @param connectance_NM Niche Model parameter specifying the expected connectivity of the network, which determines for a given node the niche space window within which it attaches to every other node. Defaults to 0.2
#' 
#' @param force_connected Binary argument determining if the newly grown node has to be connected to the existing network. Defaults to FALSE, to prevent rare computational slow-downs when it is unlikely to create a connected network. Defaults to FALSE.
#' 
#' @param directed Binary variable determining if the network is directed, resulting in off-diagonal asymmetry in the adjacency matrix. Defaults to TRUE.
#' 
#' @details This function grows, one node at a time, a mixture mechanism network. As each node is added to the growing network it can attach to existing nodes by its own node-specific mechanism. A sequence of mechanism names must be provided. Note: Currently each mechanism is assumed to have a single governing parameter.
#'
#' @return An unweighted mixture mechanism adjacency matrix.
#' 
#' @references Langendorf, R. E., & Burgess, M. G. (2020). Empirically Classifying Network Mechanisms. arXiv preprint arXiv:2012.15863.
#' 
#' @examples
#' # Mechanisms
#' sequence <- c("ER", "SW", "SW", "ER", "PA")
#' network <- stir_Mixture(sequence)
#' 
#' @export

stir_Mixture <- function(matrix, sequence, stirs = 100, p_ER = 0.5, power_PA = 2, divergence_DM = 0.8, mutation_DM = 0.1, link_DM = 0.2, connectance_NM = 0.2, rewire_SW = 0.1, force_connected = FALSE, directed = TRUE) {
    # ## Primary Directory
    # pd <- "/Users/ryan/Windows/Documents/Post UCB/Research/Relativism"
    # setwd(pd)

    # ## Libraries
    # library("tidyverse")
    # library("vegan")
    # library("igraph")

    # ## Custom Functions
    # source("stir_ER.R")
    # source("stir_PA.R")
    # source("stir_DD.R")
    # source("stir_SW.R")
    # # source("grow_CC.R")

    size <- length(sequence)

    sequence_id = seq_along(sequence)

    niches <- runif(size) %>% sort()

    sequence_stir = {}
    for (i in 1:stirs) {
        sequence_stir = c(sequence_stir, sample(x = sequence_id, size = length(sequence_id), replace = FALSE))

        for (s in seq_along(sequence_stir)) {
            x = sequence_stir[s]

            # print(matrix)
            # print(sequence[x])

            if (sequence[x] == "ER") {
                matrix = stir_ER(matrix = matrix, x = x, p = p_ER)
            } else if (sequence[x] == "PA") {
                matrix = stir_PA(matrix = matrix, x = x, power = power_PA)
            } else if (sequence[x] == "DM") {
                matrix = stir_DM(matrix = matrix, x = x, divergence = divergence_DM, mutation = mutation_DM, link = link_DM, force_connected = force_connected)
            } else if (sequence[x] == "NM") {
                matrix = stir_NM(matrix = matrix, x = x, niches = niches, connectance = connectance_NM, directed = directed)
            } else if (sequence[x] == "SW") {
                matrix = stir_SW(matrix = matrix, x = x, rewire = rewire_SW)
            } else {
                print("ERROR: Model not specified.")
            }
        }

    }

    return(matrix)
}

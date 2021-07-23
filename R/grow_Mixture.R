#' @title Grow a Mixture Mechanism Network
#'
#' @description Creates a network by iteratively adding nodes, each capable of attaching to existing nodes according to a user-specified mechanism.
#'
#' @param sequence A vector of mechanism names corresponding to the mechanisms each node acts in accordance with. Needs to be the same length as the number of nodes in the network. Note that the first two mechanisms are irrelevant because the first two nodes default to connecting to each other. Currently supported mechanisms: "ER" (Erdos-Renyi random), "PA", (Preferential Attachment), "DD", (Duplication and Divergence), "DM" (Duplication and Mutation), "SW", (Small-World), and "NM" (Niche Model).
#' 
#' @param niches Used by the Niche Model to determine which nodes interact. Needs to be a vector of the same length as the number of nodes, and range between zero and one.
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
#' @param retcon Binary variable determining if already existing nodes can attach to new nodes. Defaults to FALSE.
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
#' # Import netcom
#' library(netcom)
#' 
#' # Mechanisms
#' sequence <- c("ER", "SW", "SW", "ER", "PA")
#' network <- grow_Mixture(sequence)
#' 
#' @export

grow_Mixture <- function(sequence, niches, p_ER = 0.5, power_PA = 2, divergence_DD = 0.1, link_DD = 0, divergence_DM = 0.1, mutation_DM = 0.1, link_DM = 0, rewire_SW = 0.1, connectance_NM = 0.2, retcon = FALSE, directed = TRUE) {
    size <- length(sequence)
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

    for (x in 3:size) {
        if (sequence[x] == "ER") {
            matrix = grow_ER(matrix, x, p = p_ER, retcon = retcon, directed = TRUE) #directed)
        } else if (sequence[x] == "PA") {
            matrix = grow_PA(matrix, x, power = power_PA, retcon = retcon, directed = TRUE) #directed)
        } else if (sequence[x] == "DD") {
            matrix = grow_DD(matrix, x, divergence = divergence_DD, link = link_DD, directed = FALSE) #directed)
        } else if (sequence[x] == "DM") {
            matrix = grow_DM(matrix, x, divergence = divergence_DM, mutation = mutation_DM, link = link_DM, directed = FALSE) #directed)
        } else if (sequence[x] == "SW") {
            matrix = grow_SW(matrix, x, rewire = rewire_SW, retcon = retcon, directed = FALSE) #directed)
        } else if (sequence[x] == "NM") {
            matrix = grow_NM(matrix, x, connectance = connectance_NM, niches = niches, retcon = retcon, directed = TRUE) #directed) 
        } else {
            print("ERROR: Model not specified.")
        }
    }

    return(matrix)
}

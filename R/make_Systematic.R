#' @title Systematically Make Networks
#'
#' @description Creates a list of networks that systematically spans mechanisms and their respective parameters.
#'
#' @param net_size Number of nodes in the network.
#' 
#' @param neighborhood The range of nodes that form connected communities. Note: This implementation results in overlap of communities.
#' 
#' @param directed Whether the target network is directed. Defaults to TRUE.
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list"). Defaults to "matrix".
#' 
#' @param mechanism_kind Either "canonical" or "grow" can be used to simulate networks. If "grow" is used, note that here it will only simulate pure mixtures made of a single mechanism. Defaults to "canonical".
#' 
#' @param resolution The first step is to find the version of each process most similar to the target network. This parameter sets the number of parameter values to search across. Decrease to improve performance, but at the cost of accuracy. Defaults to 100.
#' 
#' @param resolution_min = The minimum parameter value to consider. Zero is not used because in many processes it results in degenerate systems (e.g. entirely unconnected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process. Defaults to 0.01.
#' 
#' @param resolution_max The maximum parameter value to consider. One is not used because in many processes it results in degenerate systems (e.g. entirely connected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process. Defaults to 0.99.
#' 
#' @param reps Defaults to 3. The number of networks to simulate for each parameter. More replicates increases accuracy by making the estimation of the parameter that produces networks most similar to the target network less idiosyncratic.
#' 
#' @param processes Defaults to c("ER", "PA", "DD", "SW", "NM"). Vector of process abbreviations. Currently only the default five are supported. Future versions will accept user-defined network-generating functions and associated parameters. ER = Erdos-Renyi random. PA = Preferential Attachment. DD = Duplication and Divergence. SW = Small World. NM = Niche Model.
#' 
#' @param power_max Defaults to 5. The maximum power of attachment in the Preferential Attachment process (PA).
#' 
#' @param connectance_max Defaults to 0.5. The maximum connectance parameter for the Niche Model.
#' 
#' @param divergence_max Defaults to 0.5. The maximum divergence parameter for the Duplication and Divergence/Mutation mechanisms.
#' 
#' @param mutation_max Defaults to 0.5. The maximum mutation parameter for the Duplication and Mutation mechanism.
#' 
#' @param canonical Defautls to FALSE. If TRUE the mechanisms are directed or undirected in accordance with their canonical forms. This negates the value of `directed`.
#' 
#' @param cores = Defaults to 1. The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param verbose = Defaults to TRUE. Whether to print all messages.
#'
#' @details Produces ground-truthing network data.
#'
#' @return A list. The first element contains the networks. The second contains their corresponding parameters.
#' 
#' @references Langendorf, R. E., & Burgess, M. G. (2020). Empirically Classifying Network Mechanisms. arXiv preprint arXiv:2012.15863.
#' 
#' @examples
#' make_Systematic(net_size = 10)
#' 
#' @export

make_Systematic <- function(net_size, neighborhood, directed = TRUE, net_kind = "matrix", mechanism_kind = "canonical", resolution = 100, resolution_min = 0.01, resolution_max = 0.99, reps = 3, processes = c("ER", "PA", "DM", "SW", "NM"), power_max = 5, connectance_max = 0.5, divergence_max = 0.5, mutation_max = 0.5, canonical = FALSE, cores = 1, verbose = TRUE) {

    return_list <- switch(mechanism_kind,
                            "canonical" = make_Systematic_canonical(net_size = net_size,
                                                                    net_kind = net_kind,
                                                                    resolution = resolution,
                                                                    resolution_min = resolution_min,
                                                                    resolution_max = resolution_max,
                                                                    reps = reps,
                                                                    processes = processes,
                                                                    power_max = power_max,
                                                                    connectance_max = connectance_max,
                                                                    divergence_max = divergence_max,
                                                                    mutation_max = mutation_max,
                                                                    cores = cores,
                                                                    directed = directed,
                                                                    verbose = verbose),
                            "grow" = make_Systematic_mixture(net_size = net_size,
                                                            net_kind = net_kind,
                                                            resolution = resolution,
                                                            resolution_min = resolution_min,
                                                            resolution_max = resolution_max,
                                                            reps = reps,
                                                            processes = processes,
                                                            power_max = power_max,
                                                            connectance_max = connectance_max,
                                                            divergence_max = divergence_max,
                                                            mutation_max = mutation_max,
                                                            cores = cores,
                                                            directed = directed,
                                                            verbose = verbose))

    return(return_list)
}
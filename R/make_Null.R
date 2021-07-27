#' @title Mechanism Null Distributions
#'
#' @description Creates a null distribution for a mechanism and parameter combination.
#'
#' @param input_network The network for which to create a null distribution.
#' 
#' @param net_size Number of nodes in the network.
#' 
#' @param iters Number of replicates in the null distribution. Note that length(null_dist) = ((iters^2)-iters)/2.
#' 
#' @param neighborhood The range of nodes that form connected communities. Note: This implementation results in overlap of communities.
#' 
#' @param directed Whether the target network is directed.
#' 
#' @param DD_kind = A vector of network properties to be used to compare networks.
#' 
#' @param DD_weight = A vector of weights for the relative importance of the network properties in DD_kind being used to compare networks. Should be the same length as DD_kind.
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list"). Defaults to "matrix".
#' 
#' @param mechanism_kind Either "canonical" or "grow" can be used to simulate networks. If "grow" is used, note that here it will only simulate pure mixtures made of a single mechanism.
#' 
#' @param process Name of mechanism. Currently only "ER", "PA", "DD", "DM" "SW", and "NM" are supported. Future versions will accept user-defined network-generating functions and associated parameters. ER = Erdos-Renyi random. PA = Preferential Attachment. DD = Duplication and Divergence. DM = Duplication and Mutation. SW = Small World. NM = Niche Model.
#' 
#' @param parameter Parameter in the governing mechanism.
#' 
#' @param method This determines the method used to compare networks at the heart of the classification. Currently "DD" (Degree Distribution) and "align" (the align function which compares networks by the entropy of diffusion on them) are supported. Future versions will allow user-defined methods.
#'
#' @param size_different If there is a difference in the size of the networks used in the null distribution. Defaults to FALSE.
#' 
#' @param resolution_min The minimum parameter value to consider. Zero is not used because in many processes it results in degenerate systems (e.g. entirely unconnected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process. Defaults to 0.01.
#' 
#' @param resolution_max The maximum parameter value to consider. One is not used because in many processes it results in degenerate systems (e.g. entirely connected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process. Defaults to 0.99.
#' 
#' @param power_max Defaults to 5. The maximum power of attachment in the Preferential Attachment process (PA).
#' 
#' @param connectance_max Defaults to 0.5. The maximum connectance parameter for the Niche Model.
#' 
#' @param divergence_max Defaults to 0.5. The maximum divergence parameter for the Duplication and Divergence/Mutation mechanisms.
#' 
#' @param best_fit_sd Defaults to 0.01. Standard Deviation used to simulate networks with a similar but not identical best fit parameter. This is important because simulating networks with the identical parameter artificially inflates the false negative rate by assuming the best fit parameter is the true parameter. For large resolution and reps values this will become true, but also computationally intractable for realistically large systems.
#' 
#' @param max_norm Binary variable indicating if each network property should be normalized so its max value (if a node-level property) is one. Defaults to FALSE.
#' 
#' @param cause_orientation The orientation of directed adjacency matrices. Defaults to "row".
#' 
#' @param cores Defaults to 1. The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param verbose Defaults to FALSE. Whether to print all messages.
#'
#' @details Produces ground-truthing network data.
#'
#' @return A list. The first element contains the networks. The second contains their corresponding parameters.
#' 
#' @references Langendorf, R. E., & Burgess, M. G. (2020). Empirically Classifying Network Mechanisms. arXiv preprint arXiv:2012.15863.
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' make_Systematic(net_size = 10)
#' 
#' @export

make_Null <- function(input_network, net_kind, mechanism_kind, process, parameter, net_size, iters, method, neighborhood, DD_kind, DD_weight, directed, resolution_min = 0.01, resolution_max = 0.99, power_max = 5, connectance_max = 0.5, divergence_max = 0.5, best_fit_sd = 0, cores = 1, size_different = FALSE, cause_orientation = "row", max_norm = FALSE, verbose = FALSE) {

    null_list <- switch(mechanism_kind,
                        "canonical" = make_Null_canonical(input_network = input_network,
                                                            net_kind = net_kind,
                                                            DD_kind = DD_kind,
                                                            DD_weight = DD_weight,
                                                            process = process, 
                                                            parameter = parameter,
                                                            power_max = power_max,
                                                            connectance_max = connectance_max,
                                                            divergence_max = divergence_max,
                                                            net_size = net_size, 
                                                            iters = iters, ## Note: length(null_dist) = ((iters^2)-iters)/2
                                                            method = method,
                                                            neighborhood = max(1, round(0.1 * net_size)),
                                                            best_fit_sd = best_fit_sd,
                                                            cores = cores,
                                                            directed = directed,
                                                            size_different = size_different,
                                                            cause_orientation = cause_orientation,
                                                            max_norm = max_norm,
                                                            verbose = verbose),
                        "grow" = make_Null_mixture(input_network = input_network,
                                                    net_kind = net_kind,
                                                    DD_kind = DD_kind,
                                                    DD_weight = DD_weight,
                                                    process = process, 
                                                    parameter = parameter,
                                                    power_max = power_max,
                                                    connectance_max = connectance_max,
                                                    divergence_max = divergence_max,
                                                    net_size = net_size, 
                                                    iters = iters, ## Note: length(null_dist) = ((iters^2)-iters)/2
                                                    method = method,
                                                    neighborhood = max(1, round(0.1 * net_size)),
                                                    best_fit_sd = best_fit_sd,
                                                    cores = cores,
                                                    directed = directed,
                                                    size_different = size_different,
                                                    cause_orientation = cause_orientation,
                                                    max_norm = max_norm,
                                                    verbose = verbose))

    return(null_list)
}

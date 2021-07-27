#' @title Empirical parameterization
#'
#' @description Helper function to find the best fitting version of a mechanism by searching across its parameter space
#'
#' @param parameter The parameter being tested for its ability to generate networks alike the input `network`.
#' 
#' @param process Name of mechanism. Currently only "ER", "PA", "DD", "DM" "SW", and "NM" are supported. Future versions will accept user-defined network-generating functions and associated parameters. ER = Erdos-Renyi random. PA = Preferential Attachment. DD = Duplication and Divergence. DM = Duplication and Mutation. SW = Small World. NM = Niche Model.
#' 
#' @param network The network being compared to a hypothesized `process` with a given `parameter` value.
#' 
#' @param net_size Number of nodes in the network.
#' 
#' @param directed Whether the target network is directed.
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list").
#' 
#' @param mechanism_kind Either "canonical" or "grow" can be used to simulate networks. If "grow" is used, note that here it will only simulate pure mixtures made of a single mechanism.
#' 
#' @param resolution The first step is to find the version of each process most similar to the target network. This parameter sets the number of parameter values to search across. Decrease to improve performance, but at the cost of accuracy.
#' 
#' @param resolution_min = The minimum parameter value to consider. Zero is not used because in many processes it results in degenerate systems (e.g. entirely unconnected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process.
#' 
#' @param resolution_max The maximum parameter value to consider. One is not used because in many processes it results in degenerate systems (e.g. entirely connected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process.
#' 
#' @param reps The number of networks to simulate for each parameter. More replicates increases accuracy by making the estimation of the parameter that produces networks most similar to the target network less idiosyncratic.
#' 
#' @param power_max The maximum power of attachment in the Preferential Attachment process (PA).
#' 
#' @param connectance_max The maximum connectance parameter for the Niche Model.
#' 
#' @param divergence_max The maximum divergence parameter for the Duplication and Divergence/Mutation mechanisms.
#' 
#' @param mutation_max The maximum mutation parameter for the Duplication and Mutation mechanism.
#' 
#' @param cores The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param method This determines the method used to compare networks at the heart of the classification. Currently "DD" (Degree Distribution) and "align" (the align function which compares networks by the entropy of diffusion on them) are supported. Future versions will allow user-defined methods.
#' 
#' @param cause_orientation The orientation of directed adjacency matrices.
#' 
#' @param DD_kind A vector of network properties to be used to compare networks.
#' 
#' @param DD_weight Weights of each network property in DD_kind. Defaults to 1, which is equal weighting for each property.
#' 
#' @param max_norm Binary variable indicating if each network property should be normalized so its max value (if a node-level property) is one.
#' 
#' @param best_fit_kind How to aggregate the stochastic replicates of the process + parameter combination.
#' 
#' @param verbose Defaults to TRUE. Whether to print all messages.
#' 
#' @details Note: Currently each process is assumed to have a single governing parameter.
#'
#' @return A number measuring how different the input network is from the parameter + process combination.
#' 
#' @references Langendorf, R. E., & Burgess, M. G. (2020). Empirically Classifying Network Mechanisms. arXiv preprint arXiv:2012.15863.
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' # Adjacency matrix
#' size <- 10
#' network <- matrix(sample(c(0,1), size = size^2, replace = TRUE), nrow = size, ncol = size)
#' 
#' # Calculate how similar the input network is to Small-World networks with 
#' # a rewiring probability of 0.28.
#' best_fit_optim(
#'      parameter = 0.28, 
#'      process = "SW", 
#'      network = network, 
#'      net_size = 12, 
#'      net_kind = "matrix", 
#'      mechanism_kind = "grow", 
#'      resolution = 100, 
#'      resolution_min = 0.01, 
#'      resolution_max = 0.99, 
#'      reps = 3, 
#'      power_max = 5, 
#'      connectance_max = 0.5, 
#'      divergence_max = 0.5, 
#'      mutation_max = 0.5, 
#'      cores = 1, 
#'      directed = TRUE, 
#'      method = "DD", 
#'      cause_orientation = "row", 
#'      DD_kind = c(
#'          "in", "out", "entropy_in", "entropy_out", 
#'          "clustering_coefficient", "page_rank", "communities"
#'      ), 
#'      DD_weight = 1, 
#'      max_norm = FALSE,
#'      verbose = FALSE
#' )
#' 
#' @export

best_fit_optim <- function(parameter, process, network, net_size, net_kind, mechanism_kind, resolution, resolution_min, resolution_max, reps, power_max, connectance_max, divergence_max, mutation_max, cores, directed, method, cause_orientation, DD_kind, DD_weight, max_norm, best_fit_kind = "avg", verbose = FALSE) {
    ## Create object with list of networks and table of corresponding parameters
    state_space <- make_Systematic(net_size = net_size,
                                net_kind = net_kind,
                                mechanism_kind = mechanism_kind,
                                resolution = 1,
                                resolution_min = parameter,
                                resolution_max = parameter,
                                reps = reps,
                                processes = process,
                                power_max = power_max,
                                connectance_max = connectance_max,
                                divergence_max = divergence_max,
                                mutation_max = mutation_max,
                                cores = cores,
                                directed = directed,
                                verbose = verbose)

    networks <- state_space$networks
    parameters <- state_space$parameters

    D_target <- compare_Target(target = network, 
                            networks = networks, 
                            net_kind = net_kind,
                            method = method, 
                            cause_orientation = cause_orientation, 
                            DD_kind = DD_kind, 
                            DD_weight = DD_weight,
                            max_norm = max_norm, 
                            cores = cores, 
                            verbose = verbose)

    if (best_fit_kind == "avg") {
        optim_fit <- mean(D_target)
    } else if (best_fit_kind == "min") {
        optim_fit <- min(D_target)
    } else if (best_fit_kind == "max") {
        optim_fit <- max(D_target)
    } else if (best_fit_kind == "median") {
        optim_fit <- stats::median(D_target)
    } else {
        stop("best_fit_kind must be `avg`, `median`, `min`, or `max`.")
    }

    if (verbose) { print(paste0("best_fit_optim: ", parameter, " => ", optim_fit)) }

    return(optim_fit)
}

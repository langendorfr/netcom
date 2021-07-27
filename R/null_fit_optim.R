#' @title Empirical parameterization via null distributions
#'
#' @description Helper function to find the best fitting version of a mechanism by searching across the null distributions associated with a process + parameter combination.
#'
#' @param parameter The parameter being tested for its ability to generate networks alike the input `network`.
#' 
#' @param process Name of mechanism. Currently only "ER", "PA", "DD", "DM" "SW", and "NM" are supported. Future versions will accept user-defined network-generating functions and associated parameters. ER = Erdos-Renyi random. PA = Preferential Attachment. DD = Duplication and Divergence. DM = Duplication and Mutation. SW = Small World. NM = Niche Model.
#' 
#' @param network The network being compared to a hypothesized `process` with a given `parameter` value.
#' 
#' @param net_size Number of nodes in the network.
#' 
#' @param iters Number of replicates in the null distribution. Note that length(null_dist) = ((iters^2)-iters)/2.
#' 
#' @param neighborhood The range of nodes that form connected communities. Note: This implementation results in overlap of communities.
#' 
#' @param directed Whether the target network is directed.
#' 
#' @param DD_kind A vector of network properties to be used to compare networks.
#' 
#' @param DD_weight A vector of weights for the relative importance of the network properties in DD_kind being used to compare networks. Should be the same length as DD_kind.
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list").
#' 
#' @param mechanism_kind Either "canonical" or "grow" can be used to simulate networks. If "grow" is used, note that here it will only simulate pure mixtures made of a single mechanism.
#' 
#' @param method This determines the method used to compare networks at the heart of the classification. Currently "DD" (Degree Distribution) and "align" (the align function which compares networks by the entropy of diffusion on them) are supported. Future versions will allow user-defined methods.
#'
#' @param size_different If there is a difference in the size of the networks used in the null distribution.
#' 
#' @param power_max The maximum power of attachment in the Preferential Attachment process (PA).
#' 
#' @param connectance_max The maximum connectance parameter for the Niche Model.
#' 
#' @param divergence_max The maximum divergence parameter for the Duplication and Divergence/Mutation mechanisms.
#' 
#' @param best_fit_sd Standard Deviation used to simulate networks with a similar but not identical best fit parameter. This is important because simulating networks with the identical parameter artificially inflates the false negative rate by assuming the best fit parameter is the true parameter. For large resolution and reps values this will become true, but also computationally intractable for realistically large systems.
#' 
#' @param max_norm Binary variable indicating if each network property should be normalized so its max value (if a node-level property) is one.
#' 
#' @param cause_orientation The orientation of directed adjacency matrices.
#' 
#' @param cores The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param null_dist_trim = Number between zero and one that determines how much of each network comparison distribution (unknown network compared to simulated networks, simulated networks compared to each other) should be used. Prevents p-value convergence with large sample sizes. Defaults to 1, which means all comparisons are used (no trimming).
#' 
#' @param ks_dither The KS test cannot compute exact p-values when every pairwise network distance is not unique. Adding small amounts of noise makes each distance unique. We are not aware of a study on the impacts this has on accuracy so it is set to zero by default.
#' 
#' @param ks_alternative Governs the KS test. Assuming best_fit_sd is not too large, this can be set to "greater" because the target network cannot be more alike identically simulated networks than they are to each other. In practice we have found "greater" and "less" produce numerical errors. Only "two.sided", "less", and "greater" are supported through stats::ks.test().
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
#' null_fit_optim(
#'      parameter = 0.28, 
#'      process = "SW", 
#'      network = network, 
#'      net_size = 12, 
#'      iters = 20,
#'      neighborhood = max(1, round(0.1 * net_size)),
#'      net_kind = "matrix", 
#'      mechanism_kind = "grow", 
#'      power_max = 5, 
#'      connectance_max = 0.5, 
#'      divergence_max = 0.5, 
#'      cores = 1, 
#'      directed = TRUE, 
#'      method = "DD", 
#'      size_different = FALSE,
#'      cause_orientation = "row", 
#'      DD_kind = c(
#'          "in", "out", "entropy_in", "entropy_out", 
#'          "clustering_coefficient", "page_rank", "communities"
#'      ), 
#'      DD_weight = 1, 
#'      best_fit_sd = 0,
#'      max_norm = FALSE,
#'      null_dist_trim = 0,
#'      ks_dither = 0,
#'      ks_alternative = "two.sided",
#'      verbose = FALSE
#' )
#' 
#' @export

null_fit_optim <- function(parameter, process, network, net_size, iters, neighborhood, directed, DD_kind, DD_weight, net_kind, mechanism_kind, method, size_different, power_max, connectance_max, divergence_max, best_fit_sd, max_norm, cause_orientation, cores, null_dist_trim, ks_dither, ks_alternative, verbose = FALSE) {

    best_fit <- parameter

    null_dist <- make_Null(input_network = network,
                            net_kind = net_kind,
                            mechanism_kind = mechanism_kind,
                            DD_kind = DD_kind,
                            DD_weight = DD_weight,
                            process = process, 
                            parameter = best_fit,
                            power_max = power_max,
                            connectance_max = connectance_max,
                            divergence_max = divergence_max,
                            net_size = net_size, 
                            iters = iters, ## Note: length(null_dist) = ((iters^2)-iters)/2
                            method = method,
                            neighborhood = neighborhood,
                            best_fit_sd = best_fit_sd,
                            cores = cores,
                            directed = directed,
                            size_different = size_different,
                            cause_orientation = cause_orientation,
                            max_norm = max_norm,
                            verbose = verbose)

    ## Distance matrix is symmetric so only use the lower triangular values
    null_dist_network <- null_dist$D_null[1,-1]
    null_dist_process <- null_dist$D_null[-1,-1]
    null_dist_process = null_dist_process[lower.tri(null_dist_process, diag = FALSE)]

    if ((null_dist_trim > 0) & (null_dist_trim < 1)) {
        null_dist_network = null_dist_network[1:round(length(null_dist_network)*null_dist_trim)]
        null_dist_process = null_dist_process[1:round(length(null_dist_process)*null_dist_trim)]
    }

    ### KS test
    ## alternative should be "less" because less than expected just means best_fit_sd is too big, but I was getting weird p-values of 1 so left the default as "two.sided"
    ks_test <- stats::ks.test(x = null_dist_network + stats::rnorm(n = length(null_dist_network), mean = 0, sd = ks_dither), 
                                y = null_dist_process + stats::rnorm(n = length(null_dist_process), mean = 0, sd = ks_dither), 
                                alternative = ks_alternative) %>% suppressWarnings() 

    p_value <- ks_test$p.value

    if (verbose) { print(paste0(parameter, " => ", 1 - p_value)) }

    return(1 - p_value)
}

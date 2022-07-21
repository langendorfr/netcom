#' @title Mechanistic Network Classification
#'
#' @description Tests a network against hypothetical generating processes using a comparative network inference.
#'
#' @param network The network to be classified.
#' 
#' @param directed Defaults to TRUE. Whether the target network is directed.
#' 
#' @param method This determines the method used to compare networks at the heart of the classification. Currently "DD" (Degree Distribution) and "align" (the align function which compares networks by the entropy of diffusion on them) are supported. Future versions will allow user-defined methods. Defaults to "DD".
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list"). Defaults to "matrix".
#' 
#' @param DD_kind = A vector of network properties to be used to compare networks. Defaults to "all", which is the average of the in- and out-degrees.
#' 
#' @param DD_weight = Weights of each network property in DD_kind. Defaults to 1, which is equal weighting for each property.
#' 
#' @param cause_orientation = The orientation of directed adjacency matrices. Defaults to "row".
#' 
#' @param max_norm Binary variable indicating if each network property should be normalized so its max value (if a node-level property) is one. Defaults to FALSE.
#'
#' @param resolution Defaults to 100. The first step is to find the version of each process most similar to the target network. This parameter sets the number of parameter values to search across. Decrease to improve performance, but at the cost of accuracy.
#' 
#' @param resolution_min Defaults to 0.01. The minimum parameter value to consider. Zero is not used because in many processes it results in degenerate systems (e.g. entirely unconnected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process.
#' 
#' @param resolution_max Defaults to 0.99. The maximum parameter value to consider. One is not used because in many processes it results in degenerate systems (e.g. entirely connected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process.
#' 
#' @param connectance_max = Defaults to 0.5. The maximum connectance parameter for the Niche Model.
#' 
#' @param divergence_max = Defaults to 0.5. The maximum divergence parameter for the Duplication and Divergence/Mutation mechanisms.
#' 
#' @param mutation_max = Defaults to 0.5. The maximum mutation parameter for the Duplication and Mutation mechanism.
#' 
#' @param reps Defaults to 3. The number of networks to simulate for each parameter. More replicates increases accuracy by making the estimation of the parameter that produces networks most similar to the target network less idiosyncratic.
#' 
#' @param processes Defaults to c("ER", "PA", "DD", "SW", "NM"). Vector of process abbreviations. Currently only the default five are supported. Future versions will accept user-defined network-generating functions and associated parameters. ER = Erdos-Renyi random. PA = Preferential Attachment. DD = Duplication and Divergence. SW = Small World. NM = Niche Model.
#' 
#' @param power_max Defaults to 5. The maximum power of attachment in the Preferential Attachment process (PA).
#' 
#' @param null_reps Defaults to 50. The number of best fit networks to simulate that will be used to create a null distribution of distances between networks within the given process, which will then be used to test if the target network appears unusually distant from them and therefore likely not governed by that process.
#' 
#' @param best_fit_kind Defaults to "avg". If null_reps is more than 1, the fit of each parameter has to be an aggregate statistic of the fit of all the null_reps networks. Must be `avg`, `median`, `min`, or `max`.
#' 
#' @param best_fit_sd Defaults to 0.01. Standard Deviation used to simulate networks with a similar but not identical best fit parameter. This is important because simulating networks with the identical parameter artificially inflates the false negative rate by assuming the best fit parameter is the true parameter. For large resolution and reps values this will become true, but also computationally intractable for realistically large systems.
#' 
#' @param ks_dither Defaults to 0. The KS test cannot compute exact p-values when every pairwise network distance is not unique. Adding small amounts of noise makes each distance unique. We are not aware of a study on the impacts this has on accuracy so it is set to zero by default.
#' 
#' @param ks_alternative Defaults to "two.sided". Governs the KS test. Assuming best_fit_sd is not too large, this can be set to "greater" because the target network cannot be more alike identically simulated networks than they are to each other. In practice we have found "greater" and "less" produce numerical errors. Only "two.sided", "less", and "greater" are supported through stats::ks.test().
#' 
#' @param cores Defaults to 1. The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param size_different = If there is a difference in the size of the networks used in the null distribution. Defaults to FALSE.
#' 
#' @param null_dist_trim = Number between zero and one that determines how much of each network comparison distribution (unknown network compared to simulated networks, simulated networks compared to each other) should be used. Prevents p-value convergence with large sample sizes. Defaults to 1, which means all comparisons are used (no trimming).
#' 
#' @param verbose Defaults to TRUE. Whether to print all messages.
#' 
#' @details Note: Currently each process is assumed to have a single governing parameter.
#'
#' @return A dataframe with 3 columns and as many rows as processes being tested (5 by default). The first column lists the processes. The second lists the p-value on the null hypothesis that the target network did come from that row's process. The third column gives the estimated parameter for that particular process.
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
#' # Classify this network
#' # This can take several minutes to run
#' \donttest{classify(network, processes = c("ER", "PA", "DM", "SW", "NM"))}
#' 
#' @export
#' 
#' @importFrom rlang .data

classify_Systematic <- function(network, directed = FALSE, method = "DD", net_kind = "matrix", DD_kind = c("in", "out", "entropy_in", "entropy_out", "clustering_coefficient", "page_rank", "communities", "motifs_3", "motifs_4", "eq_in", "eq_out", "eq_entropy_in", "eq_entropy_out", "eq_clustering_coefficient", "eq_page_rank", "eq_communities", "eq_motifs_3", "eq_motifs_4"), DD_weight = c(0.0735367966, 0.0739940162, 0.0714523761, 0.0708156931, 0.0601296752, 0.0448072016, 0.0249793608, 0.0733125084, 0.0697029389, 0.0504358835, 0.0004016029, 0.0563752664, 0.0561878218, 0.0540490099, 0.0504347104, 0.0558106667, 0.0568270319, 0.0567474398), cause_orientation = "row", max_norm = FALSE, resolution = 100, resolution_min = 0.01, resolution_max = 0.99, reps = 3, processes = c("ER", "PA", "DM", "SW", "NM"), power_max = 5, connectance_max = 0.5, divergence_max = 0.5, mutation_max = 0.5, null_reps = 50, best_fit_kind = "avg", best_fit_sd = 1e-2, ks_dither = 0, ks_alternative = "two.sided", cores = 1, size_different = FALSE, null_dist_trim = 1, verbose = TRUE) {

    ## Matrix input checks
    if (net_kind == "matrix") {
        ## Check that the network is square
        if (nrow(network) != ncol(network)) {
            stop("Input network must be a square matrix.")
        }

        ## If there are row and column names, check that they are ordered identically
        ## If there are none, assume they are
        if((!is.null(rownames(network))) & (!is.null(colnames(network)))) {
            if (!identical(rownames(network), colnames(network))) {
                stop("Row & Column names must have the same order.")
            }
        }
    }

    ## Network size
    if (net_kind == "matrix") {
        net_size <- nrow(network)
    } else if (net_kind == "list") {
        net_size <- c(network[,1], network[,2]) %>% unique() %>% length()
    } else {
        stop("Unknown network kind. Must be `list` or `matrix`.")
    }

    ## Create object with list of networks and table of corresponding parameters
    state_space <- make_Systematic(net_size = net_size,
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
                                   verbose = verbose)

    networks <- state_space$networks
    parameters <- state_space$parameters

    D_target <- compare_Target(target = network, 
                               networks = networks, 
                            #    net_size = net_size,
                               net_kind = net_kind,
                               method = method, 
                               cause_orientation = cause_orientation, 
                               DD_kind = DD_kind, 
                               DD_weight = DD_weight,
                               max_norm = max_norm, 
                            #    size_different = size_different,
                               cores = cores, 
                               verbose = verbose)
    
    # print("D_target")
    # print(D_target)

    parameters_scored <- dplyr::mutate(parameters, Distance = D_target)

    ## Check each process one at a time
    pvalues <- rep(NA, length(processes))
    p_estimates <- rep(NA, length(processes))
    for (p in seq_along(processes)) {
        if (verbose) {print(paste0("Checking if the network is ", processes[p], "."))}

        parameters_scored_process <- dplyr::filter(parameters_scored, .data$Process == processes[p])

        ## Use the min of the average
        ## Even for a given process and parameter there are many possible networks

        if (best_fit_kind == "avg") {
            parameters_scored_process = parameters_scored_process %>% dplyr::group_by(.data$Process, .data$Parameter_Value) %>% dplyr::summarize(Distance = mean(.data$Distance), .groups = "drop")
        } else if (best_fit_kind == "min") {
            parameters_scored_process = parameters_scored_process %>% dplyr::group_by(.data$Process, .data$Parameter_Value) %>% dplyr::summarize(Distance = min(.data$Distance), .groups = "drop")
        } else if (best_fit_kind == "max") {
            parameters_scored_process = parameters_scored_process %>% dplyr::group_by(.data$Process, .data$Parameter_Value) %>% dplyr::summarize(Distance = max(.data$Distance), .groups = "drop")
        } else if (best_fit_kind == "median") {
            parameters_scored_process = parameters_scored_process %>% dplyr::group_by(.data$Process, .data$Parameter_Value) %>% dplyr::summarize(Distance = stats::median(.data$Distance), .groups = "drop")
        } else {
            stop("best_fit_kind must be `avg`, `min`, or `max`.")
        }

        best_fit <- dplyr::filter(parameters_scored_process, .data$Distance == min(parameters_scored_process$Distance))

        ## Assuming there are multiple best fits, pick one randomly
        best_fit = best_fit[sample(1:nrow(best_fit), size = 1), ]

        # print("best_fit")
        # print(best_fit)

        null_dist <- make_Null(input_network = network,
                               net_kind = net_kind,
                               DD_kind = DD_kind,
                               DD_weight = DD_weight,
                               process = best_fit$Process, 
                               parameter = best_fit$Parameter_Value,
                               power_max = power_max,
                               connectance_max = connectance_max,
                               divergence_max = divergence_max,
                               net_size = net_size, 
                               iters = null_reps, ## Note: length(null_dist) = ((iters^2)-iters)/2
                               method = method,
                               neighborhood = max(1, round(0.1 * net_size)),
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

        ## alternative should be "less" because less than expected just means best_fit_sd is too big, but I was getting weird p-values of 1 so left the default as "two.sided"
        ks_test <- stats::ks.test(x = null_dist_network + stats::rnorm(n = length(null_dist_network), mean = 0, sd = ks_dither), 
                                  y = null_dist_process + stats::rnorm(n = length(null_dist_process), mean = 0, sd = ks_dither), 
                                  alternative = ks_alternative) %>% suppressWarnings() 

        # p_value <- 1 - sum(best_fit$Distance > null_dist) / length(null_dist)
        p_value <- ks_test$p.value

        pvalues[p] = p_value
        # p_estimates[p] = stats::weighted.mean(x = null_dist$parameters, w = 1/(null_dist$D_null[1,-1]))
        p_estimates[p] = stats::weighted.mean(x = null_dist$parameters, w = exp(-null_dist$D_null[1,-1]) )
    }

    return_tbl <- tibble::tibble(process = processes, 
                         p_value = pvalues, #round(pvalues, 5),
                         par_estimate = p_estimates) #round(p_estimates, 5))
    

    # ## Warn about strangely close best_fit distances
    # for (row in 1:nrow(return_tbl)) {
    #     if (return_tbl$p_value[row] > 0.95) {
    #         warning(paste0("Network appears strangely similar to process ", return_tbl$Process[row]))
    #     }
    # }


    return(return_tbl)
}
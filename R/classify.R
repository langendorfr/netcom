#' @title Mechanistic Network Classification
#'
#' @description Tests a network against hypothetical generating processes using a comparative network inference.
#'
#' @param network The network to be classified.
#' 
#' @param directed Whether the target network is directed. If missing this will be inferred by the symmetry of the input network.
#' 
#' @param method This determines the method used to compare networks at the heart of the classification. Currently "DD" (Degree Distribution) and "align" (the align function which compares networks by the entropy of diffusion on them) are supported. Future versions will allow user-defined methods. Defaults to "DD".
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list"). Defaults to "matrix".
#' 
#' @param mechanism_kind Either "canonical" or "grow" can be used to simulate networks. If "grow" is used, note that here it will only simulate pure mixtures made of a single mechanism. Defaults to "canonical".
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
#' @param test Defaults to "empirical". The test used to distinguish the null distribution of comparisons between the network being classified and the networks simulated according to a hypothesized mechanism(s), with a particular best-fitting parameter. "empirical" finds how many simulated networks were on average farther from each other than the network being classified is. "KS" uses a KS test. "WMWU" uses a Wilcoxon-Mann-Whitney-U test.
#' 
#' @param best_fit_finder Defaults to "systematic". Determines how the best-fitting parameter of each mechanism specified in processes is found. "systematic" tries every parameter value from resolution_min to resolution_max with a step size of resolution_max - resolution_min / resolution. "optim_L-BFGS-B" uses the L-BFGS-B optimizer in the optimx package. "optim_GenSA" uses the GenSA optimizer in the GenSA package.
#' 
#' @param best_fit_kind Defaults to "avg". If null_reps is more than 1, the fit of each parameter has to be an aggregate statistic of the fit of all the null_reps networks. Must be `avg`, `median`, `min`, or `max`.
#' 
#' @param best_fit_sd Defaults to 0. Standard Deviation used to simulate networks with a similar but not identical best fit parameter. This is important because simulating networks with the identical parameter can artificially inflate the false negative rate by assuming the best fit parameter is the true parameter. For large resolution and reps values this will become true, but can be computationally intractable for realistically large systems.
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
#' @param verbose Defaults to FALSE. Whether to print all messages.
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
#' \dontrun{classify(network, processes = c("ER", "PA", "DM", "SW", "NM"))}
#' 
#' @export

classify <- function(network, directed, method = "DD", net_kind = "matrix", mechanism_kind = "canonical", DD_kind = c("in", "out", "entropy_in", "entropy_out", "clustering_coefficient", "page_rank", "communities", "motifs_3", "motifs_4", "eq_in", "eq_out", "eq_entropy_in", "eq_entropy_out", "eq_clustering_coefficient", "eq_page_rank", "eq_communities", "eq_motifs_3", "eq_motifs_4"), DD_weight = c(0.0735367966, 0.0739940162, 0.0714523761, 0.0708156931, 0.0601296752, 0.0448072016, 0.0249793608, 0.0733125084, 0.0697029389, 0.0504358835, 0.0004016029, 0.0563752664, 0.0561878218, 0.0540490099, 0.0504347104, 0.0558106667, 0.0568270319, 0.0567474398), cause_orientation = "row", max_norm = FALSE, resolution = 100, resolution_min = 0.01, resolution_max = 0.99, reps = 3, processes = c("ER", "PA", "DM", "SW", "NM"), test = "empirical", best_fit_finder = "systematic", power_max = 5, connectance_max = 0.5, divergence_max = 0.5, mutation_max = 0.5, null_reps = 50, best_fit_kind = "avg", best_fit_sd = 0, ks_dither = 0, ks_alternative = "two.sided", cores = 1, size_different = FALSE, null_dist_trim = 1, verbose = FALSE) {

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

    ## Determine if the network is directed
    directed_inferred <- !isSymmetric(network)
    if (missing(directed)) {
        directed <- directed_inferred
    }

    ## Warn if the network appears to be un/directed differently from what was specified by the user
    if (directed_inferred != directed) {
        if (directed) {
            warning_directed <- "is "
        } else {
            warning_directed <- "is not "
        }

        warning(paste0("You specified that the network ", warning_directed, "directed, but the symmetry of the network implies this is not the case." ))
    }

    ## Network size
    if (net_kind == "matrix") {
        net_size <- nrow(network)
    } else if (net_kind == "list") {
        net_size <- c(network[,1], network[,2]) %>% unique() %>% length()
    } else {
        stop("Unknown network kind. Must be `list` or `matrix`.")
    }

    ## Check each process one at a time
    pvalues <- rep(NA, length(processes))
    p_estimates <- rep(NA, length(processes))
    for (p in seq_along(processes)) {
        if (verbose) {print(paste0("Checking if the network is ", processes[p], "."))}

        if (best_fit_finder == "systematic") {
            possible_pars <- seq(from = resolution_min, to = resolution_max, length.out = resolution)
            possible_pars_fits <- rep(NA, length(possible_pars))
            for (par in seq_along(possible_pars)) {
                if (verbose) {print(paste0("Comparing the network to ", processes[p], " with parameter ", possible_pars[par], "."))}
                possible_pars_fits[par] = best_fit_optim(parameter = possible_pars[par], process = processes[p], best_fit_kind = best_fit_kind, network = network, net_size = net_size, net_kind = net_kind, mechanism_kind = mechanism_kind, resolution = resolution, resolution_min = resolution_min, resolution_max = resolution_max, reps = reps, power_max = power_max, connectance_max = connectance_max, divergence_max = divergence_max, mutation_max = mutation_max, cores = cores, directed = directed, method = method, cause_orientation = cause_orientation, DD_kind = DD_kind, DD_weight = DD_weight, max_norm = max_norm)
            }

            ## Use the first best fitting parameter in case there are multiple. Do not aggregate these in case they are different optima
            best_fit <- possible_pars[which(possible_pars_fits == min(possible_pars_fits))][1]

        } else if (best_fit_finder == "optim_L-BFGS-B") {
            optim_fit <- optimx::optimx(par = mean(c(resolution_min, resolution_max)), ## Start in the midde of the range of possible parameter values
                                        fn = best_fit_optim, 
                                        method = "L-BFGS-B", 
                                        lower = resolution_min, 
                                        upper = resolution_max) #[1,1]
            best_fit = optim_fit[1,1]

        } else if (best_fit_finder == "optim_GenSA") {
            optim_fit <- GenSA::GenSA(fn = best_fit_optim,  #null_fit_optim, #best_fit_optim, 
                                        lower = resolution_min, 
                                        upper = resolution_max,
                                        control = list(max.call = 1e2, maxit = 1e2, smooth = FALSE))
            best_fit = optim_fit$par
        } else {
            stop("best_fit_finder must be 'systematic', 'optim_L-BFGS-B', or 'optim_GenSA'.")
        }

        null_dist <- make_Null(input_network = network,
                               net_kind = net_kind,
                               mechanism_kind = mechanism_kind,
                               DD_kind = DD_kind,
                               DD_weight = DD_weight,
                               process = processes[p], 
                               parameter = best_fit,
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

        ## Distances from target network to hypothesized mechanism with the best fit parameter
        null_dist_network <- null_dist$D_null[1,-1]

        if (test == "empirical") {
            ## Best fitting network simulated from the mechanism with the best fit parameter
            best_fit_null_network <- which(null_dist_network == min(null_dist_network))[1]
            ## +1 because left out the first row for the target network in calculating `null_dist_network`
            null_dist_process <- null_dist$D_null[1+best_fit_null_network, -c(1, 1+best_fit_null_network)]


            ## Percentile of average distance to networks simulated from the mechanism with the best fit parameter
            avg_null_distances <- (rowSums(null_dist$D_null)/nrow(null_dist$D_null))
            p_value <- sum(avg_null_distances[-1] > avg_null_distances[1])/(length(avg_null_distances)-1)

                                        avg_null_distances <- apply(null_dist$D_null, 2, mean)
                                        var_null_distances <- apply(null_dist$D_null, 2, stats::var)
                                        avg_null_distances = avg_null_distances * var_null_distances
                                        p_value <- sum(avg_null_distances[-1] > avg_null_distances[1])/(length(avg_null_distances)-1)
                                        if (p_value > 0.5) {
                                            p_value = 1 - p_value
                                        }

        } else if (test == "KS") {
            # Distance matrix is symmetric so only use the lower triangular values
            null_dist_process <- null_dist$D_null[-1,-1]
            null_dist_process = null_dist_process[lower.tri(null_dist_process, diag = FALSE)]

            if ((null_dist_trim > 0) & (null_dist_trim < 1)) {
                null_dist_network = null_dist_network[1:round(length(null_dist_network)*null_dist_trim)]
                null_dist_process = null_dist_process[1:round(length(null_dist_process)*null_dist_trim)]
            }

            ## KS test
            # alternative should be "less" because less than expected just means best_fit_sd is too big, but I was getting weird p-values of 1 so left the default as "two.sided"
            ks_test <- stats::ks.test(x = null_dist_network + stats::rnorm(n = length(null_dist_network), mean = 0, sd = ks_dither), 
                                    y = null_dist_process + stats::rnorm(n = length(null_dist_process), mean = 0, sd = ks_dither), 
                                    alternative = ks_alternative) %>% suppressWarnings() 

            # p_value <- 1 - sum(best_fit$Distance > null_dist) / length(null_dist)
            p_value <- ks_test$p.value

        } else if (test == "WMWU") {
            # Distance matrix is symmetric so only use the lower triangular values
            null_dist_process <- null_dist$D_null[-1,-1]
            null_dist_process = null_dist_process[lower.tri(null_dist_process, diag = FALSE)]

            if ((null_dist_trim > 0) & (null_dist_trim < 1)) {
                null_dist_network = null_dist_network[1:round(length(null_dist_network)*null_dist_trim)]
                null_dist_process = null_dist_process[1:round(length(null_dist_process)*null_dist_trim)]
            }
        
            ## WMWU test
            p_value <- stats::wilcox.test(null_dist_network, null_dist_process, alternative = "greater")$p.value

        } else {
            stop("test must be 'empirical', 'KS', or 'WMWU'.")
        }

        pvalues[p] = p_value
        p_estimates[p] = stats::weighted.mean(x = null_dist$parameters, w = exp(-null_dist$D_null[1,-1]) )
    }

    return_tbl <- tibble::tibble(process = processes, 
                         par_estimate = p_estimates,
                         p_value = pvalues)

    return(return_tbl)
}
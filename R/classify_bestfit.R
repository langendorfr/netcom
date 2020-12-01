#' @title 
#'
#' @description 
#'
#' @param 
#'
#' @details 
#'
#' @return 
#' 
#' @references 
#' 
#' @examples
#' 
#' @export

classify_bestfit <- function(network, method = "DD", net_kind = "list", resolution = 100, resolution_min = 0.01, resolution_max = 0.99, reps = 3, processes = c("ER", "PA", "DD", "SW", "NM"), power_max = 5, null_reps = 50, best_fit_sd = 1e-2, ks_dither = 0, ks_alternative = "two.sided", cores = 1, directed = TRUE, verbose = TRUE) {

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
                                   cores = cores,
                                   directed = directed,
                                   verbose = verbose)

    networks <- state_space$networks
    parameters <- state_space$parameters

    D_target <- compare_Target(target = network, 
                               networks = networks, 
                               net_size = net_size,
                               net_kind = net_kind,
                               method = method, 
                               cause_orientation = "row", 
                               DD_kind = "out", 
                               max_norm = FALSE, 
                               cores = 1, 
                               verbose = FALSE)

    parameters_scored <- dplyr::mutate(parameters, Distance = D_target)

    ## Check each process one at a time
    pvalues <- rep(NA, length(processes))
    p_estimates <- rep(NA, length(processes))
    bestfit_values <- rep(NA, length(processes))
    for (p in seq_along(processes)) {
        if (verbose == TRUE) {
            print(paste0("Checking if the network is ", processes[p], "."))
        }

        parameters_scored_process <- dplyr::filter(parameters_scored, Process == processes[p])

        ## Use the min of the average
        ## Even for a given process and parameter there are many possible networks
        parameters_scored_process = parameters_scored_process %>% group_by(Process, Parameter_Value) %>% summarize(Distance = mean(Distance), .groups = "drop")
        best_fit <- dplyr::filter(parameters_scored_process, Distance == min(parameters_scored_process$Distance))

        ## Assuming there are multiple best fits, pick one randomly
        best_fit = best_fit[sample(1:nrow(best_fit), size = 1), ]

        bestfit_values[p] = best_fit
#######################
        pvalues[p] = p_value
        # p_estimates[p] = stats::weighted.mean(x = null_dist$parameters, w = 1/(null_dist$D_null[1,-1]))
        p_estimates[p] = stats::weighted.mean(x = null_dist$parameters, w = exp(-null_dist$D_null[1,-1]) )
    }

    return_tbl <- tibble(process = processes, 
                        #  p_value = pvalues, #round(pvalues, 5),
                        #  par_estimate = p_estimates) #round(p_estimates, 5)
                        bestfit = bestfit_values
                         )
    

    # ## Warn about strangely close best_fit distances
    # for (row in 1:nrow(return_tbl)) {
    #     if (return_tbl$p_value[row] > 0.95) {
    #         warning(paste0("Network appears strangely similar to process ", return_tbl$Process[row]))
    #     }
    # }


    return(return_tbl)
}
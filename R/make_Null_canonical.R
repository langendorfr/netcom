#' @title Mechanism Null Distributions
#'
#' @description Creates a null distribution for a mechanism and parameter combination.
#'
#' @param input_network The network for which to create a null distribution.
#' 
#' @param net_size Number of nodes in the network.
#' 
#' @param neighborhood The range of nodes that form connected communities. Note: This implementation results in overlap of communities.
#' 
#' @param directed Whether the target network is directed. Defaults to TRUE.
#' 
#' @param DD_kind = A vector of network properties to be used to compare networks. Defaults to "all", which is the average of the in- and out-degrees.
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list"). Defaults to "matrix".
#' 
#' @param resolution The first step is to find the version of each process most similar to the target network. This parameter sets the number of parameter values to search across. Decrease to improve performance, but at the cost of accuracy. Defaults to 100.
#' 
#' @param resolution_min = The minimum parameter value to consider. Zero is not used because in many processes it results in degenerate systems (e.g. entirely unconnected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process. Defaults to 0.01.
#' 
#' @param resolution_max The maximum parameter value to consider. One is not used because in many processes it results in degenerate systems (e.g. entirely connected networks). Currently process agnostic. Future versions will accept a vector of values, one for each process. Defaults to 0.99.
#' 
#' @param reps Defaults to 3. The number of networks to simulate for each parameter. More replicates increases accuracy by making the estimation of the parameter that produces networks most similar to the target network less idiosyncratic.
#' 
#' @param process Name of mechanism. Currently only "ER", "PA", "DD", "DM" "SW", and "NM" are supported. Future versions will accept user-defined network-generating functions and associated parameters. ER = Erdos-Renyi random. PA = Preferential Attachment. DD = Duplication and Divergence. DM = Duplication and Mutation. SW = Small World. NM = Niche Model.
#' 
#' @param parameter Parameter in the governing mechanism.
#' 
#' @param method This determines the method used to compare networks at the heart of the classification. Currently "DD" (Degree Distribution) and "align" (the align function which compares networks by the entropy of diffusion on them) are supported. Future versions will allow user-defined methods.
#'
#' @param size_different If there is a difference in the size of the networks used in the null distribution. Defaults to FALSE.
#' 
#' @param DD_resize = If networks being compared are a different size, this parameter determines if upscaling "larger" or downscaling "smaller" occurs. Unlikely to be relevant here. Defaults to "smaller".
#' 
#' @param power_max = Defaults to 5. The maximum power of attachment in the Preferential Attachment process (PA).
#' 
#' @param connectance_max = Defaults to 0.5. The maximum connectance parameter for the Niche Model.
#' 
#' @param divergence_max = Defaults to 0.5. The maximum divergence parameter for the Duplication and Divergence/Mutation mechanisms.
#' 
#' @param mutation_max = Defaults to 0.5. The maximum mutation parameter for the Duplication and Mutation mechanism.
#' 
#' @param best_fit_sd Defaults to 0.01. Standard Deviation used to simulate networks with a similar but not identical best fit parameter. This is important because simulating networks with the identical parameter artificially inflates the false negative rate by assuming the best fit parameter is the true parameter. For large resolution and reps values this will become true, but also computationally intractable for realistically large systems.
#' 
#' @param max_norm Binary variable indicating if each network property should be normalized so its max value (if a node-level property) is one. Defaults to FALSE.
#' 
#' @param cause_orientation = The orientation of directed adjacency matrices. Defaults to "row".
#' 
#' @param cores = Defaults to 1. The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param verbose = Defaults to FALSE. Whether to print all messages.
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

make_Null <- function(input_network, net_kind, process, parameter, net_size, iters, method, neighborhood, DD_kind, DD_weight, resolution_min = 0.01, resolution_max = 0.99, directed = TRUE, power_max = 5, connectance_max = 0.5, divergence_max = 0.5, best_fit_sd = 0, cores = 1, size_different = FALSE, cause_orientation = "row", max_norm = FALSE, DD_resize = "smaller", verbose = FALSE) {
    ## Primary Directory
    # pd <- "/Users/ryan/Windows/Documents/Post UCB/Research/Relativism"
    # setwd(pd)

    # ## Libraries
    # library("tidyverse")
    # library("vegan")
    # library("igraph")

    # ## Custom Functions
    # source("grow_ER.R")
    # source("grow_PA.R")
    # source("grow_DD.R")
    # source("grow_SW.R")
    # source("grow_CM.R")

    if (!(net_kind %in% c("matrix", "list"))) {
        stop("Unknown net_kind. Must be `list` or `matrix`.")
    }

    ## Build list of networks with same process and parameter (+ some noise depending on the best_fit_sd parameter)
    networks <- list()
    parameters <- {}
    
    ## Start with the network being classified
    networks[[1]] = input_network

    ## +1 because include the network being classified
    for (counter in 2:(iters+1)) {



                ## Add network to growing list, made from same process and parameter
                if (process == "ER") {
                    # directed = TRUE

                    p_ER <- parameter + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    p_ER = min(p_ER, resolution_max)
                    p_ER = max(p_ER, resolution_min)
                    parameters = c(parameters, p_ER)
                    
                    net <- igraph::sample_gnp(n = net_size, 
                                              p = p_ER, 
                                              directed = directed, 
                                              loops = FALSE)

                    if (net_kind == "matrix") {
                        mat <- igraph::as_adj(net, 
                                              type = "both", 
                                              edges = TRUE, 
                                              names = TRUE,
                                              sparse = FALSE)
                        ## igraph puts the edge id in the matrix element
                        mat[which(mat != 0)] = 1
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        # edgelist <- igraph::as_edgelist(net)
                        # edgelist <- mat %>% as.matrix() %>% reshape2::melt() %>% dplyr::filter(value != 0)
                        edgelist <- net %>% igraph::as.directed(mode = "mutual") %>% igraph::as_edgelist(names = TRUE)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                } else if (process == "PA") {
                    # directed = TRUE

                    power_PA <- (parameter * power_max) + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    power_PA = min(power_PA, resolution_max * power_max)
                    power_PA = max(power_PA, resolution_min * power_max)
                    parameters = c(parameters, power_PA)

                    net <- igraph::sample_pa(n = net_size, 
                                             power = power_PA,
                                             directed = directed,
                                             m = 1, #NULL, 
                                             out.dist = NULL, 
                                             out.seq = NULL, 
                                             out.pref = FALSE, 
                                             zero.appeal = 1,
                                             algorithm = "psumtree",
                                             start.graph = NULL)

                    if (net_kind == "matrix") {
                        mat <- igraph::as_adj(net, 
                                            type = "both", 
                                            edges = TRUE, 
                                            names = TRUE,
                                            sparse = FALSE)
                        ## igraph puts the edge id in the matrix element
                        mat[which(mat != 0)] = 1
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- net %>% igraph::as.directed(mode = "mutual") %>% igraph::as_edgelist(names = TRUE)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                } else if (process == "DD") {
                    # directed = FALSE

                    divergence_DD <- (parameter * divergence_max) + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    divergence_DD = min(divergence_DD, resolution_max)
                    divergence_DD = max(divergence_DD, resolution_min)
                    parameters = c(parameters, divergence_DD)

                    net <- make_DD(size = net_size, 
                                   net_kind = net_kind, 
                                   divergence = divergence_DD, 
                                   directed = directed)

                    networks[[counter]] = net

                } else if (process == "DM") {
                    # directed = FALSE

                    divergence_DM <- (parameter  * divergence_max) + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    divergence_DM = min(divergence_DM, resolution_max)
                    divergence_DM = max(divergence_DM, resolution_min)
                    parameters = c(parameters, divergence_DM)

                    mutation_DM <- divergence_DM

                    net <- make_DM(size = net_size, 
                                   net_kind = net_kind, 
                                   divergence = divergence_DM,
                                   mutation = mutation_DM, 
                                   directed = directed)

                    networks[[counter]] = net

                } else if (process == "SW") {
                    # directed = FALSE

                    rewire_SW <- parameter + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    rewire_SW = min(rewire_SW, resolution_max)
                    rewire_SW = max(rewire_SW, resolution_min)
                    parameters = c(parameters, rewire_SW)

                    ## SW neighborhood parameter based on net_size if missing
                    if (missing(neighborhood)) {
                        neighborhood = max(1, round(0.1 * net_size))
                    }

                    net <- make_SW(size = net_size, 
                                   net_kind = net_kind, 
                                   rewire = rewire_SW, 
                                   neighborhood = neighborhood, 
                                   directed = directed)

                    networks[[counter]] <- net

                } else if (process == "NM") {
                    # directed = TRUE

                    connectance_NM <- (parameter * connectance_max) + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    connectance_NM = min(connectance_NM, resolution_max * connectance_max)
                    connectance_NM = max(connectance_NM, resolution_min * connectance_max)
                    parameters = c(parameters, connectance_NM)

                    niches <- runif(net_size) # %>% sort()
                    net <- make_NM(size = net_size,
                                   net_kind = net_kind,
                                   niches = niches, 
                                   connectance = connectance_NM, 
                                   directed = directed,
                                   grow = TRUE)
                    networks[[counter]] <- net

                } else {
                    stop("An unknown process was included in the simulation.")
                }








    } ## for (counter in 1:iters) {

    ## Compare uuid networks to get distribution of their distances
    # if (verbose == TRUE) { print("Mapping network state space.") }

    D_null <- compare(networks = networks, 
                      method = method, 
                      net_kind = net_kind,
                    #   net_size = net_size,
                      cause_orientation = cause_orientation,
                      DD_kind = DD_kind,
                      DD_weight = DD_weight,
                      DD_resize = DD_resize,
                      max_norm = max_norm,
                      size_different = size_different,
                      cores = cores, 
                      verbose = verbose)
    
    # D_null  %>% gplots::heatmap.2(trace = "none", Rowv = FALSE, Colv = FALSE, dendrogram = "none", col = colorRampPalette(c("white", "black"))(n = 299))
    
    null_list <- list(D_null = D_null,
                      parameters = parameters)

    return(null_list)
}

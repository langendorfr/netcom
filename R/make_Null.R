#' @title 
#'
#' @description 
#'
#' @param 
#'
#' @details Only makes directed networks
#'
#' @return 
#' 
#' @references 
#' 
#' @examples
#' 
#' @export

make_Null <- function(input_network, net_kind, process, parameter, net_size, iters, method, best_fit_sd = 0, cores = 1, directed = TRUE, verbose = TRUE) {
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

    ## Build list of networks with same process and parameter (+ some noise depending on the best_fit_sd parameter)
    networks <- list()
    parameters <- {}
    
    ## Start with the network being classified
    networks[[1]] = input_network

    ## +1 because include the network being classified
    for (counter in 2:(iters+1)) {



                ## Add network to growing list, made from same process and parameter
                if (process == "ER") {
                    p_ER <- parameter + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    p_ER = min(p_ER, resolution_max)
                    p_ER = max(p_ER, resolution_min)
                    parameters = c(parameters, p_ER)
                    
                    net <- igraph::sample_gnp(n = net_size, 
                                            p = p_ER, 
                                            directed = directed, 
                                            loops = FALSE)

                    mat <- igraph::as_adj(net, 
                                        type = "both", 
                                        edges = TRUE, 
                                        names = TRUE,
                                        sparse = FALSE)
                    ## igraph puts the edge id in the matrix element
                    mat[which(mat != 0)] = 1

                    if (net_kind == "matrix") {
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- mat %>% as.matrix() %>% reshape2::melt() %>% dplyr::filter(value != 0)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                } else if (process == "PA") {
                    power_PA <- parameter + (rnorm(n = 1, mean = 0, sd = best_fit_sd) * power_max)
                    power_PA = min(power_PA, resolution_max * power_max)
                    power_PA = max(power_PA, resolution_min * power_max)
                    parameters = c(parameters, power_PA)

                    net <- igraph::sample_pa(n = net_size, 
                                            power = power_PA, 
                                            m = 1, #NULL, 
                                            out.dist = NULL, 
                                            out.seq = NULL, 
                                            out.pref = FALSE, 
                                            zero.appeal = 1, 
                                            directed = directed, 
                                            algorithm = "psumtree",
                                            start.graph = NULL)

                    mat <- igraph::as_adj(net, 
                                        type = "both", 
                                        edges = TRUE, 
                                        names = TRUE,
                                        sparse = FALSE)
                    ## igraph puts the edge id in the matrix element
                    mat[which(mat != 0)] = 1

                    if (net_kind == "matrix") {
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- mat %>% as.matrix() %>% reshape2::melt() %>% dplyr::filter(value != 0)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                } else if (process == "DD") {
                    divergence_DD <- parameter + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    divergence_DD = min(divergence_DD, resolution_max)
                    divergence_DD = max(divergence_DD, resolution_min)
                    parameters = c(parameters, divergence_DD)                    
                    
                    net <- nx$duplication_divergence_graph(n=net_size, 
                                                        p=divergence_DD)
                    
                    if (net_kind == "matrix") {
                        mat <- nx$to_numpy_matrix(G=net, 
                                                weight='weight', 
                                                nonedge=0.0)
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- nx$to_pandas_edgelist(G=net)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                } else if (process == "SW") {
                    rewire_SW <- parameter + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    rewire_SW = min(rewire_SW, resolution_max)
                    rewire_SW = max(rewire_SW, resolution_min)
                    parameters = c(parameters, rewire_SW)

                    net <- igraph::sample_smallworld(dim = 1, 
                                                    size = net_size, 
                                                    nei = max(2, round(0.1 * net_size)), 
                                                    p = rewire_SW, 
                                                    loops = FALSE, 
                                                    multiple = FALSE)

                    mat <- igraph::as_adj(net, 
                                        type = "both", 
                                        edges = TRUE, 
                                        names = TRUE,
                                        sparse = FALSE)
                    ## igraph puts the edge id in the matrix element
                    mat[which(mat != 0)] = 1

                    if (net_kind == "matrix") {
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- mat %>% as.matrix() %>% reshape2::melt() %>% dplyr::filter(value != 0)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                } else if (process == "NM") {
                    connectance_NM <- parameter + rnorm(n = 1, mean = 0, sd = best_fit_sd)
                    connectance_NM = min(connectance_NM, resolution_max)
                    connectance_NM = max(connectance_NM, resolution_min)
                    parameters = c(parameters, connectance_NM)

                    niches <- runif(net_size) # %>% sort()
                    net <- make_NM(size = net_size,
                                   net_kind = net_kind,
                                   niches = niches, 
                                   connectance = connectance_NM, 
                                   directed = directed)
                    networks[[counter]] <- net

                } else {
                    stop("An unknown process was included in the simulation.")
                }








    } ## for (counter in 1:iters) {

    ## Compare uuid networks to get distribution of their distances
    # if (verbose == TRUE) { print("Mapping network state space.") }

    D_null <- compare(networks = networks, 
                      method = method, 
                      net_size = net_size,
                      cause_orientation = "row", 
                      DD_kind = "out", 
                      DD_resize = "smaller", 
                      max_norm = FALSE, 
                      cores = cores, 
                      verbose = FALSE)

    null_list <- list(D_null = D_null,
                      parameters = parameters)

    return(null_list)
}

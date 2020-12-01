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

make_Null <- function(input_network, net_kind, process, parameter, net_size, iters, method, neighborhood, directed, DD_kind, power_max = 5, connectance_max = 0.5, divergence_max = 0.5, best_fit_sd = 0, cores = 1, size_different = FALSE, cause_orientation = "row", max_norm = FALSE, DD_resize = "smaller", verbose = FALSE) {
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
                    directed = TRUE

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
                    directed = TRUE

                    power_PA <- parameter + (rnorm(n = 1, mean = 0, sd = best_fit_sd) * power_max)
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
                    directed = FALSE

                    divergence_DD <- parameter + (rnorm(n = 1, mean = 0, sd = best_fit_sd) * divergence_max)
                    divergence_DD = min(divergence_DD, resolution_max)
                    divergence_DD = max(divergence_DD, resolution_min)
                    parameters = c(parameters, divergence_DD)

                    net <- make_DD(size = net_size, 
                                   net_kind = net_kind, 
                                   divergence = divergence_DD, 
                                   directed = directed)

                    networks[[counter]] = net

                } else if (process == "DM") {
                    directed = FALSE

                    divergence_DM <- parameter + (rnorm(n = 1, mean = 0, sd = best_fit_sd) * divergence_max)
                    divergence_DM = min(divergence_DM, resolution_max)
                    divergence_DM = max(divergence_DM, resolution_min)
                    parameters = c(parameters, divergence_DM)

                    mutation_DM <- divergence_DM

                    net <- make_DM(size = net_size, 
                                   net_kind = net_kind, 
                                   divergence = divergence_DM,
                                   mutation = 0, #mutation_DM, 
                                   directed = directed)

                    networks[[counter]] = net

                } else if (process == "SW") {
                    directed = FALSE

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
                    directed = TRUE

                    connectance_NM <- parameter + (rnorm(n = 1, mean = 0, sd = best_fit_sd) * connectance_max)
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

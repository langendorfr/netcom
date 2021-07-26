#' @title Compare Networks Many-to-Many
#'
#' @description Compares one network to a list of many networks.
#' 
#' @param networks The networks being compared to the target network
#' 
#' @param net_kind If the network is an adjacency matrix ("matrix") or an edge list ("list"). Defaults to "matrix".
#' 
#' @param method This determines the method used to compare networks at the heart of the classification. Currently "DD" (Degree Distribution) and "align" (the align function which compares networks by the entropy of diffusion on them) are supported. Future versions will allow user-defined methods. Defaults to "DD".
#' 
#' @param cause_orientation = The orientation of directed adjacency matrices. Defaults to "row".
#' 
#' @param DD_kind = A vector of network properties to be used to compare networks. Defaults to "all", which is the average of the in- and out-degrees.
#' 
#' @param DD_weight = Weights of each network property in DD_kind. Defaults to 1, which is equal weighting for each property.
#' 
#' @param max_norm Binary variable indicating if each network property should be normalized so its max value (if a node-level property) is one. Defaults to FALSE.
#' 
#' @param size_different Defaults to FALSE. If TRUE, will ensure the node-level properties being compared are vectors of the same length, which is accomplished using splines.
#'
#' @param cores Defaults to 1. The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param diffusion_sampling Base of the power to use to nonlinearly sample the diffusion kernels if method = "align". Defaults to 2.
#' 
#' @param diffusion_limit Number of markov steps in the diffusion kernels if method = "align". Defaults to 10.
#' 
#' @param verbose Defaults to TRUE. Whether to print all messages.
#' 
#' @details Note: Currently each process is assumed to have a single governing parameter.
#'
#' @return A square matrix with dimensions equal to the number of networks being compared, where the ij element is the comparison of networks i and j.
#' 
#' @references Langendorf, R. E., & Burgess, M. G. (2020). Empirically Classifying Network Mechanisms. arXiv preprint arXiv:2012.15863.
#' 
#' @examples
#' # Import netcom
#' library(netcom)
#' 
#' # Adjacency matrix
#' size <- 10
#' comparisons <- 50
#' networks <- list()
#' for (net in 1:comparisons) {
#'      networks[[net]] = matrix(
#'          sample(
#'              c(0,1), 
#'              size = size^2, 
#'              replace = TRUE), 
#'          nrow = size,
#'          ncol = size)
#' }
#' compare(networks = networks)
#' 
#' @export
#' 
#' @importFrom parallel makeCluster
#' 
#' @importFrom doParallel registerDoParallel

compare <- function(networks, net_kind = "matrix", method = "DD", cause_orientation = "row", DD_kind = "all", DD_weight = 1, max_norm = FALSE, size_different = FALSE, cores = 1, diffusion_sampling = 2, diffusion_limit = 10, verbose = FALSE) {

    if ((length(DD_weight) == 1) & (length(DD_kind) > 1)) {
        DD_weight = rep(DD_weight, length(DD_kind)) / length(DD_kind)
    }

    ## Network alignment
    if (method == "align") {
        ## Pairwise compare input network with each network in state_space
        if (cores == 1) {
            D_netcom <- matrix(NA, 
                        nrow = length(networks), 
                        ncol = length(networks))
            for (net_1 in seq_along(networks)) {
                for (net_2 in seq_along(networks)) {
                    if(net_2 > net_1) {
                        # if (verbose == TRUE) { print(c(net_1, net_2)) }

                        alignment <- align(
                            networks[[net_1]], 
                            networks[[net_2]], 
                            base = diffusion_sampling,
                            max_duration = diffusion_limit,
                            characterization = "entropy", 
                            normalization = FALSE
                        )

                        D_netcom[net_1, net_2] = alignment$score
                        D_netcom[net_2, net_1] = alignment$score
                    }
                }
            }
            diag(D_netcom) = 0

        } else { ## Cores > 1
            ## Use parallelization
            cluster <- parallel::makeCluster(cores, outfile = "")
            doParallel::registerDoParallel(cluster)

            ## Define %dopar% relative to its package so can use normally
            `%dopar%` <- foreach::`%dopar%`

            ## Network Alignment
            D_netcom <- foreach::foreach (net_1 = 1:length(networks), .combine = rbind, .packages = c("tibble", "dplyr", "netcom")) %dopar% {
                j_output <- rep(NA, length(networks))
                for (net_2 in 1:length(networks)) {
                    #  if (verbose == TRUE) { print(c(net_1, net_2)) }
                
                    j_output[net_2] <- align(
                        networks[[net_1]], 
                        networks[[net_2]], 
                        base = diffusion_sampling,
                        max_duration = diffusion_limit, 
                        characterization = "entropy", 
                        normalization = FALSE)$score
                }
                
                return(j_output) 
            }

            parallel::stopCluster(cluster)
        }

        return_matrix <- D_netcom

###################################################################################################################

    ## Degree Distribution comparisons, using Euclidean distance
    } else if (method == "DD") {
        if (net_kind == "matrix") {
            ## NOTE: assumes row -> col orientation
            ## Do nothing if (cause_orientation == "row")

            if (cause_orientation == "col") {
                networks = lapply(networks, function(x){t(x)})
            }

            D_DD <- matrix(NA, 
                        nrow = length(networks),
                        ncol = length(networks))

            # ## Assume all networks are the same size
            # net_size <- nrow(networks[[1]])

            net_DD_list <- list()
            for (net in seq_along(networks)) {
                # if (verbose == TRUE) { print(net) }

                ## Assume all networks are the same size
                net_size <- nrow(networks[[net]])
                
                DD_combined <- list()
            
                igraph_graph <- igraph::graph_from_adjacency_matrix(networks[[net]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)
                
                equilibrium_net <- networks[[net]]
                zero_rows <- which(rowSums(equilibrium_net) == 0)
                diag(equilibrium_net)[zero_rows] = 1
                equilibrium_net = sweep(equilibrium_net, 1, Matrix::rowSums(equilibrium_net, na.rm = TRUE), FUN = "/")
                equilibrium_net = expm::`%^%`(equilibrium_net, 5)

                eq_igraph_graph <- igraph::graph_from_adjacency_matrix(equilibrium_net, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)

                if ("communities" %in% DD_kind) {
                    community_sizes <- igraph::graph_from_adjacency_matrix(networks[[net]], mode = "undirected", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA) %>% igraph::cluster_fast_greedy(weights = NULL) %>% igraph::sizes() %>% sort(decreasing = TRUE)
                }

                if ("eq_communities" %in% DD_kind) {
                    eq_community_sizes <- igraph::graph_from_adjacency_matrix(equilibrium_net, mode = "undirected", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA) %>% igraph::cluster_fast_greedy(weights = NULL) %>% igraph::sizes() %>% sort(decreasing = TRUE)
                }

                for (DD_kind_name in seq_along(DD_kind)) {
                    DD_net <- switch(
                        DD_kind[DD_kind_name],
                        "out" = networks[[net]] %>% rowSums() %>% as.numeric(),
                        "eq_out" = equilibrium_net %>% rowSums() %>% as.numeric(),
                        
                        "in" = networks[[net]] %>% colSums() %>% as.numeric(),
                        "eq_in" = equilibrium_net %>% colSums() %>% as.numeric(),
                        
                        "undirected" = as.numeric(rowSums(networks[[net]])) + as.numeric(colSums(networks[[net]])),
                        "eq_undirected" = as.numeric(rowSums(equilibrium_net)) + as.numeric(colSums(equilibrium_net)),
                        
                        "all" = 0.5 * (rbind( rowSums(networks[[net]]) + colSums(networks[[net]]) )),
                        "eq_all" = 0.5 * (rbind( rowSums(equilibrium_net) + colSums(equilibrium_net) )),
                        
                        "entropy_out" = networks[[net]] %>% t() %>% vegan::diversity() %>% as.numeric(),
                        "eq_entropy_out" = equilibrium_net %>% t() %>% vegan::diversity() %>% as.numeric(),
                        
                        "entropy_in" = networks[[net]] %>% vegan::diversity() %>% as.numeric(),
                        "eq_entropy_in" = equilibrium_net %>% vegan::diversity() %>% as.numeric(),
                        
                        "entropy_all" = c(networks[[net]] %>% t() %>% vegan::diversity() %>% as.numeric(), networks[[net]] %>% vegan::diversity() %>% as.numeric()),
                        "eq_entropy_all" = c(equilibrium_net %>% t() %>% vegan::diversity() %>% as.numeric(), equilibrium_net %>% vegan::diversity() %>% as.numeric()),
                        
                        "clustering_coefficient" = igraph::transitivity(igraph_graph, type = "weighted"),
                        "eq_clustering_coefficient" = igraph::transitivity(eq_igraph_graph, type = "weighted"),

                        "cohesion" = igraph::vertex_connectivity(igraph_graph),
                        "eq_cohesion" = igraph::vertex_connectivity(eq_igraph_graph),

                        "spectral_decomposition" = igraph::embed_adjacency_matrix(igraph_graph, no = nrow(networks[[net]])-1)$X %>% rowSums(),
                        "eq_spectral_decomposition" = igraph::embed_adjacency_matrix(eq_igraph_graph, no = nrow(equilibrium_net)-1)$X %>% rowSums(),

                        "alpha" = igraph::alpha_centrality(igraph_graph, nodes = igraph::V(igraph_graph), alpha = 10, loops = TRUE, exo = 1, weights = NULL, tol = 1e-07, sparse = FALSE),
                        "eq_alpha" = igraph::alpha_centrality(eq_igraph_graph, nodes = igraph::V(eq_igraph_graph), alpha = 10, loops = TRUE, exo = 1, weights = NULL, tol = 1e-07, sparse = FALSE),

                        "page_rank" = igraph::page_rank(igraph_graph, algo = c("prpack", "arpack", "power"), vids = igraph::V(igraph_graph), directed = TRUE, damping = 0.85, personalized = NULL, weights = NULL, options = NULL)$vector,
                        "eq_page_rank" = igraph::page_rank(eq_igraph_graph, algo = c("prpack", "arpack", "power"), vids = igraph::V(eq_igraph_graph), directed = TRUE, damping = 0.85, personalized = NULL, weights = NULL, options = NULL)$vector,

                        "communities" = rep(seq_along(community_sizes), community_sizes),
                        "eq_communities" = rep(seq_along(eq_community_sizes), eq_community_sizes),

                        "motifs_3" = igraph::motifs(igraph_graph, size = 3)[-c(1,2,4)],
                        "eq_motifs_3" = igraph::motifs(eq_igraph_graph, size = 3)[-c(1,2,4)],

                        "motifs_4" = igraph::motifs(igraph_graph, size = 4)[-c(1, 2, 3, 5, 6, 7, 10, 11, 12, 16, 23, 24, 28, 29, 34, 35, 40, 63, 121)],
                        "eq_motifs_4" = igraph::motifs(eq_igraph_graph, size = 4)[-c(1, 2, 3, 5, 6, 7, 10, 11, 12, 16, 23, 24, 28, 29, 34, 35, 40, 63, 121)],

                        stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                    )

                    ## Some network properties return NA or NaN or Inf on disconnected networks
                    if (any(is.na(DD_net)) | any(is.infinite(DD_net))) {
                        DD_net[which(is.na(DD_net))] = 0
                        DD_net[which(is.infinite(DD_net))] = 0
                    }

                    DD_net = DD_net %>% sort()
                    if (max_norm == TRUE) {
                        DD_net = DD_net / max(DD_net, na.rm = TRUE)
                    }

                    ## Some network properties return NA or NaN or Inf on disconnected networks
                    if (any(is.na(DD_net)) | any(is.infinite(DD_net))) {
                        DD_net[which(is.na(DD_net))] = 0
                        DD_net[which(is.infinite(DD_net))] = 0
                    }

                    DD_combined[[DD_kind_name]] = DD_net

                }
                
                net_DD_list[[net]] = DD_combined
            }



            for (net_1 in seq_along(networks)) {
                for (net_2 in seq_along(networks)) {
                    if(net_2 > net_1) {
                        # if (verbose == TRUE) { print(c(net_1, net_2)) }

                        # DD_difference_each <- {}
                        DD_difference_each <- rep(NA, length(DD_kind))
                        for (DD_kind_name in seq_along(DD_kind)) {
                            DD_1 <- net_DD_list[[net_1]][[DD_kind_name]]
                            DD_2 <- net_DD_list[[net_2]][[DD_kind_name]]




                            if (length(DD_1) != length(DD_2)) {
                                larger <- which(c(length(DD_1), length(DD_2)) == max(c(length(DD_1), length(DD_2))))

                                if (larger == 1) {
                                    DD_larger <- DD_1
                                    DD_smaller <- DD_2
                                } else if (larger == 2) {
                                    DD_larger <- DD_2
                                    DD_smaller <- DD_1
                                } else {
                                    stop("Error in compare(). DD_* should be either 1 or 2.")
                                }

                                Position_c <- seq(from = min(DD_larger), to = max(DD_larger), length = length(DD_smaller))
                                DD_larger = stats::spline(DD_larger, xout = Position_c)$y
                            } else {

                                ## Equal size networks so does not matter which is the smaller/larger one
                                DD_smaller <- DD_1
                                DD_larger <- DD_2
                            }



                            DD_1_vs_2 <- ((DD_smaller - DD_larger)^2) #%>% sum() %>% sqrt()

                            if (max(DD_1_vs_2, na.rm = TRUE) != 0) {
                                DD_1_vs_2 = (sum(DD_1_vs_2) / max(DD_1_vs_2, na.rm = TRUE)) * DD_weight[DD_kind_name]
                            } else {
                                DD_1_vs_2 = 0
                            }

                            DD_difference_each[DD_kind_name] = DD_1_vs_2
                        }

                        DD_difference <- sqrt(sum(DD_difference_each, na.rm = TRUE))

                        ## Add to the growing D_DD matrix
                        D_DD[net_1, net_2] = DD_difference
                        D_DD[net_2, net_1] = DD_difference
                    }
                }
            }
            diag(D_DD) = 0

            return_matrix <- D_DD

        } else if (net_kind == "list") {

            D_DD <- matrix(NA, 
                        nrow = length(networks),
                        ncol = length(networks))

            for (net_1 in seq_along(networks)) {
                for (net_2 in seq_along(networks)) {
                    if(net_2 > net_1) {
                        # if (verbose == TRUE) { print(c(net_1, net_2)) }

                        ## First network's Degree Distribution
                        comparison_net <- networks[[net_1]]

                        ## Some networks from make_Systematic() will have zero edges
                        if (nrow(comparison_net) == 0) {
                            DD_1 <- rep(0, net_size)
                        } else {
                            nodes <- c(comparison_net[,1], comparison_net[,2]) %>% unique()

                            ## Give weights of one if missing
                            if (ncol(comparison_net) == 2) {
                                comparison_net = cbind(comparison_net, 1)
                            }

                            ## NOTE: This block assumes the comparison networks are of the same size as the target network
                            DD_1 <- switch(
                                DD_kind,
                                "out" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges <- comparison_net[,1] == n
                                                DD_node <- comparison_net[relevant_edges, 3] %>% sum()
                                                DD_vector = c(DD_vector, DD_node)
                                            }
                                            DD_vector
                                        },
                                "in" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges <- comparison_net[,2] == n
                                                DD_node <- comparison_net[relevant_edges, 3] %>% sum()
                                                DD_vector = c(DD_vector, DD_node)
                                            }
                                            DD_vector
                                        },
                                "undirected" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges <- (comparison_net[,1] == n) | (comparison_net[,2] == n)
                                                intersection_eges <- (comparison_net[,1] == n) & (comparison_net[,2] == n)
                                                DD_node <- sum(comparison_net[relevant_edges, 3]) - sum(comparison_net[intersection_eges, 3])
                                                DD_vector = c(DD_vector, DD_node)
                                            }
                                            DD_vector
                                        },
                                "all" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges_out <- comparison_net[,1] == n
                                                DD_node_out <- comparison_net[relevant_edges_out, 3] %>% sum()

                                                relevent_edges_in <- comparison_net[,2] == n
                                                DD_node_in <- comparison_net[relevent_edges_in, 3] %>% sum()

                                                DD_vector = c(DD_vector, DD_node_out, DD_node_in)
                                            }
                                            DD_vector
                                        },                                        
                                stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                            )

                            DD_1 = DD_1 %>% sort()
                            if (max_norm == TRUE) {
                                DD_1 = DD_1 / max(DD_1)
                            }
                        }

                        if (size_different) {
                            Position_1 <- seq(from = 0, to = 1, length = length(DD_1))
                            Points_1 <- tibble::tibble(Position = Position_1, DD = DD_1)

                            Position_c_1 <- seq(from = 0, to = 1, length = net_size)
                            Points_c_1 <- apply(Points_1, 2, function(u) stats::spline(Position_1, u, xout = Position_c_1)$y) %>% tibble::as_tibble()
                        }




                        ## Second network's Degree Distribution
                        comparison_net <- networks[[net_2]]

                        ## Some networks from make_Systematic() will have zero edges
                        if (nrow(comparison_net) == 0) {
                            DD_2 <- rep(0, net_size)
                        } else {
                            nodes <- c(comparison_net[,1], comparison_net[,2]) %>% unique()

                            ## Give weights of one if missing
                            if (ncol(comparison_net) == 2) {
                                comparison_net = cbind(comparison_net, 1)
                            }

                            ## NOTE: This block assumes the comparison networks are of the same size as the target network
                            DD_2 <- switch(
                                DD_kind,
                                "out" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges <- comparison_net[,1] == n
                                                DD_node <- comparison_net[relevant_edges, 3] %>% sum()
                                                DD_vector = c(DD_vector, DD_node)
                                            }
                                            DD_vector
                                        },
                                "in" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges <- comparison_net[,2] == n
                                                DD_node <- comparison_net[relevant_edges, 3] %>% sum()
                                                DD_vector = c(DD_vector, DD_node)
                                            }
                                            DD_vector
                                        },
                                "undirected" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges <- (comparison_net[,1] == n) | (comparison_net[,2] == n)
                                                intersection_eges <- (comparison_net[,1] == n) & (comparison_net[,2] == n)
                                                DD_node <- sum(comparison_net[relevant_edges, 3]) - sum(comparison_net[intersection_eges, 3])
                                                DD_vector = c(DD_vector, DD_node)
                                            }
                                            DD_vector
                                        },
                                "all" = {
                                            DD_vector <- {}
                                            for (n in 1:net_size) {
                                                relevant_edges_out <- comparison_net[,1] == n
                                                DD_node_out <- comparison_net[relevant_edges_out, 3] %>% sum()

                                                relevent_edges_in <- comparison_net[,2] == n
                                                DD_node_in <- comparison_net[relevent_edges_in, 3] %>% sum()

                                                DD_vector = c(DD_vector, DD_node_out, DD_node_in)
                                            }
                                            DD_vector
                                        },
                                stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                            )

                            DD_2 = DD_2 %>% sort()
                            if (max_norm == TRUE) {
                                DD_2 = DD_2 / max(DD_2)
                            }
                        }

                        if (size_different) {
                            Position_2 <- seq(from = 0, to = 1, length = length(DD_2))
                            Points_2 <- tibble::tibble(Position = Position_2, DD = DD_2)

                            Position_c_2 <- seq(from = 0, to = 1, length = net_size)
                            Points_c_2 <- apply(Points_2, 2, function(u) stats::spline(Position_2, u, xout = Position_c_2)$y) %>% tibble::as_tibble()


                            ## Comapre the Degree Distributions of the two networks
                            DD_1 <- Points_c_1$DD
                            DD_2 <- Points_c_2$DD
                        }


                        DD_difference <- sqrt( sum( ((DD_1 - DD_2)^2) ) )

                        ## Add to the growing D_DD matrix
                        D_DD[net_1, net_2] = DD_difference
                        D_DD[net_2, net_1] = DD_difference
                    }
                }
            }
            diag(D_DD) = 0

            return_matrix <- D_DD


        
        
        
        
        
        
        
        
        
        
        
        } else {
            stop("Unknown network kind. Must be `list` or `matrix`.")
        }



    } else {
        stop("Method not supported.")
    }

    return(return_matrix)
}

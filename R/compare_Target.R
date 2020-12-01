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

compare_Target <- function(target, networks, net_size, net_kind, method, cause_orientation = "row", DD_kind = "all", DD_weight = 1, max_norm = FALSE, cores = 1, verbose = FALSE) {

    if (length(c(nrow(target), sapply(networks, nrow)) %>% unique()) == 1) {
        size_different = FALSE
    } else {
        size_different = TRUE
    }

    if ((length(DD_weight) == 1) & (length(DD_kind) > 1)) {
        DD_weight = rep(DD_weight, length(DD_kind)) / length(DD_kind)
    }

    ## Network alignment
    if (method == "align") {
        ## Pairwise compare input network with each network in state_space
        if (cores == 1) {
            if (verbose == TRUE) {
                print("Aligning networks on a single core.")
            }

            D_netcom <- rep(NA, times = length(networks))
            for (net in seq_along(networks)) {
                if (verbose == TRUE) {
                    print(paste0("Aligning ", net, " of ", length(networks), " networks."))
                }

                alignment <- netcom::align(target,
                                            networks[[net]],
                                            base = 2,
                                            characterization = "entropy",
                                            normalization = FALSE)$score
                D_netcom[net] = alignment
            }

        ## cores > 1
        } else {
            if (verbose == TRUE) {
                print(paste0("Aligning networks on ", cores, " cores."))
            }

            ## Backend
            ## outfile = "" prints to screen instead of a file
            cluster <- parallel::makeCluster(cores, outfile = "")
            doParallel::registerDoParallel(cluster)

            ## Parallelized network alignment
            D_netcom <- foreach (net = 1:length(networks), .combine = c, .packages = c("tibble", "dplyr", "netcom")) %dopar% {        
                if (verbose == TRUE) {
                    print(paste0("Aligning ", net, " of ", length(networks), " networks."))
                }

                output <- netcom::align(target,
                                        networks[[net]],
                                        base = 2,
                                        characterization = "entropy",
                                        normalization = FALSE)$score
                return(output) 
            }

            ## Stop backend
            stopCluster(cluster)
        }

        return_vector <- D_netcom

###################################################################################################################

    ## Degree Distribution comparisons, using Euclidean distance
    } else if (method == "DD") {

        if (net_kind == "matrix") {
            ## NOTE: assumes row -> col orientation
            ## Do nothing if (cause_orientation == "row")
            if (cause_orientation == "col") {
                target == t(target)
                networks = lapply(networks, function(x){t(x)})
            }

            DD_target_combined <- list()

            igraph_graph <- igraph::graph_from_adjacency_matrix(target, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)

            equilibrium_net <- target
            zero_rows <- which(rowSums(equilibrium_net) == 0)
            diag(equilibrium_net)[zero_rows] = 1
            equilibrium_net = sweep(equilibrium_net, 1, Matrix::rowSums(equilibrium_net, na.rm = TRUE), FUN = "/")
            equilibrium_net = equilibrium_net %^% 5

            eq_igraph_graph <- igraph::graph_from_adjacency_matrix(equilibrium_net, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)

            if ("communities" %in% DD_kind) {
                community_sizes <- igraph::graph_from_adjacency_matrix(target, mode = "undirected", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA) %>% igraph::cluster_fast_greedy(weights = NULL) %>% sizes() %>% sort(decreasing = TRUE)
            }

            if ("eq_communities" %in% DD_kind) {
                eq_community_sizes <- igraph::graph_from_adjacency_matrix(equilibrium_net, mode = "undirected", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA) %>% igraph::cluster_fast_greedy(weights = NULL) %>% sizes() %>% sort(decreasing = TRUE)
            }

            for (DD_kind_name in seq_along(DD_kind)) {
                DD_target <- switch(
                    DD_kind[DD_kind_name],

                    "out" = target %>% rowSums() %>% as.numeric(),
                    "eq_out" = equilibrium_net %>% rowSums() %>% as.numeric(),
                    
                    "in" = target %>% colSums() %>% as.numeric(),
                    "eq_in" = equilibrium_net %>% colSums() %>% as.numeric(),
                    
                    "undirected" = as.numeric(rowSums(target)) + as.numeric(colSums(target)),
                    "eq_undirected" = as.numeric(rowSums(equilibrium_net)) + as.numeric(colSums(equilibrium_net)),
                    
                    "all" = 0.5 * (rbind( rowSums(target) + colSums(target) )),
                    "eq_all" = 0.5 * (rbind( rowSums(equilibrium_net) + colSums(equilibrium_net) )),
                    
                    "entropy_out" = target %>% t() %>% vegan::diversity() %>% as.numeric(),
                    "eq_entropy_out" = equilibrium_net %>% t() %>% vegan::diversity() %>% as.numeric(),
                    
                    "entropy_in" = target %>% vegan::diversity() %>% as.numeric(),
                    "eq_entropy_in" = equilibrium_net %>% vegan::diversity() %>% as.numeric(),
                    
                    "entropy_all" = c(target %>% t() %>% vegan::diversity() %>% as.numeric(), target %>% vegan::diversity() %>% as.numeric()),
                    "eq_entropy_all" = c(equilibrium_net %>% t() %>% vegan::diversity() %>% as.numeric(), equilibrium_net %>% vegan::diversity() %>% as.numeric()),
                    
                    "clustering_coefficient" = igraph::transitivity(igraph_graph, type = "weighted"),
                    "eq_clustering_coefficient" = igraph::transitivity(eq_igraph_graph, type = "weighted"),

                    "cohesion" = igraph::vertex_connectivity(igraph_graph),
                    "eq_cohesion" = igraph::vertex_connectivity(eq_igraph_graph),

                    "spectral_decomposition" = igraph::embed_adjacency_matrix(igraph_graph, no = nrow(target)-1)$X %>% rowSums(),
                    "eq_spectral_decomposition" = igraph::embed_adjacency_matrix(eq_igraph_graph, no = nrow(equilibrium_net)-1)$X %>% rowSums(),

                    "alpha" = igraph::alpha_centrality(igraph_graph, nodes = V(igraph_graph), alpha = 10, loops = TRUE, exo = 1, weights = NULL, tol = 1e-07, sparse = FALSE),
                    "eq_alpha" = igraph::alpha_centrality(eq_igraph_graph, nodes = V(eq_igraph_graph), alpha = 10, loops = TRUE, exo = 1, weights = NULL, tol = 1e-07, sparse = FALSE),

                    "page_rank" = igraph::page_rank(igraph_graph, algo = c("prpack", "arpack", "power"), vids = V(igraph_graph), directed = TRUE, damping = 0.85, personalized = NULL, weights = NULL, options = NULL)$vector,
                    "eq_page_rank" = igraph::page_rank(eq_igraph_graph, algo = c("prpack", "arpack", "power"), vids = V(eq_igraph_graph), directed = TRUE, damping = 0.85, personalized = NULL, weights = NULL, options = NULL)$vector,

                    "communities" = rep(seq_along(community_sizes), community_sizes),
                    "eq_communities" = rep(seq_along(eq_community_sizes), eq_community_sizes),

                    "motifs_3" = igraph::motifs(igraph_graph, size = 3)[-c(1,2,4)],
                    "eq_motifs_3" = igraph::motifs(eq_igraph_graph, size = 3)[-c(1,2,4)],

                    "motifs_4" = igraph::motifs(igraph_graph, size = 4)[-c(1, 2, 3, 5, 6, 7, 10, 11, 12, 16, 23, 24, 28, 29, 34, 35, 40, 63, 121)],
                    "eq_motifs_4" = igraph::motifs(eq_igraph_graph, size = 4)[-c(1, 2, 3, 5, 6, 7, 10, 11, 12, 16, 23, 24, 28, 29, 34, 35, 40, 63, 121)],

                    stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                )

                ## Some network properties return NA or NaN or Inf on disconnected networks
                if (any(is.na(DD_target)) | any(is.infinite(DD_target))) {
                    DD_target[which(is.na(DD_target))] = 0
                    DD_target[which(is.infinite(DD_target))] = 0
                }
                DD_target = DD_target %>% sort()
                if (max_norm == TRUE) {
                    DD_target = DD_target / max(DD_target)
                }

                ## max(DD_target) can be zero
                if (any(is.na(DD_target)) | any(is.infinite(DD_target))) {
                    DD_target[which(is.na(DD_target))] = 0
                    DD_target[which(is.infinite(DD_target))] = 0
                }

                # DD_target_combined = c(DD_target_combined, DD_target)
                DD_target_combined[[DD_kind_name]] = DD_target
            }

            ## Just use the last DD_target, because this is the nubmer of points to use in the spline to make length(DD) = length(DD_target)
            Length <- length(DD_target)

            if (verbose == TRUE) {
                print(paste0("Comparing the degree distributions of ", length(networks), " networks."))
            }

            D_DD <- rep(NA, times = length(networks))
            for (net in seq_along(networks)) {
                if (verbose) { print(net) }
                
                DD_combined <- list()
                
                igraph_graph <- igraph::graph_from_adjacency_matrix(networks[[net]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)
                
                equilibrium_net <- networks[[net]]
                zero_rows <- which(rowSums(equilibrium_net) == 0)
                diag(equilibrium_net)[zero_rows] = 1
                equilibrium_net = sweep(equilibrium_net, 1, Matrix::rowSums(equilibrium_net, na.rm = TRUE), FUN = "/")
                equilibrium_net = equilibrium_net %^% 5

                eq_igraph_graph <- igraph::graph_from_adjacency_matrix(equilibrium_net, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)

                if ("communities" %in% DD_kind) {
                    community_sizes <- igraph::graph_from_adjacency_matrix(networks[[net]], mode = "undirected", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA) %>% igraph::cluster_fast_greedy(weights = NULL) %>% sizes() %>% sort(decreasing = TRUE)
                }

                if ("eq_communities" %in% DD_kind) {
                    eq_community_sizes <- igraph::graph_from_adjacency_matrix(equilibrium_net, mode = "undirected", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA) %>% igraph::cluster_fast_greedy(weights = NULL) %>% sizes() %>% sort(decreasing = TRUE)
                }

                for (DD_kind_name in seq_along(DD_kind)) {
                    DD <- switch(
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

                        "alpha" = igraph::alpha_centrality(igraph_graph, nodes = V(igraph_graph), alpha = 10, loops = TRUE, exo = 1, weights = NULL, tol = 1e-07, sparse = FALSE),
                        "eq_alpha" = igraph::alpha_centrality(eq_igraph_graph, nodes = V(eq_igraph_graph), alpha = 10, loops = TRUE, exo = 1, weights = NULL, tol = 1e-07, sparse = FALSE),
                        
                        "page_rank" = igraph::page_rank(igraph_graph, algo = c("prpack", "arpack", "power"), vids = V(igraph_graph), directed = TRUE, damping = 0.85, personalized = NULL, weights = NULL, options = NULL)$vector,
                        "eq_page_rank" = igraph::page_rank(eq_igraph_graph, algo = c("prpack", "arpack", "power"), vids = V(eq_igraph_graph), directed = TRUE, damping = 0.85, personalized = NULL, weights = NULL, options = NULL)$vector,

                        "communities" = rep(seq_along(community_sizes), community_sizes),
                        "eq_communities" = rep(seq_along(eq_community_sizes), eq_community_sizes),

                        "motifs_3" = igraph::motifs(igraph_graph, size = 3)[-c(1,2,4)],
                        "eq_motifs_3" = igraph::motifs(eq_igraph_graph, size = 3)[-c(1,2,4)],

                        "motifs_4" = igraph::motifs(igraph_graph, size = 4)[-c(1, 2, 3, 5, 6, 7, 10, 11, 12, 16, 23, 24, 28, 29, 34, 35, 40, 63, 121)],
                        "eq_motifs_4" = igraph::motifs(eq_igraph_graph, size = 4)[-c(1, 2, 3, 5, 6, 7, 10, 11, 12, 16, 23, 24, 28, 29, 34, 35, 40, 63, 121)],

                        stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                    )

                    ## Some network properties return NA or NaN or Inf on disconnected networks
                    if (any(is.na(DD)) | any(is.infinite(DD))) {
                        DD[which(is.na(DD))] = 0
                        DD[which(is.infinite(DD))] = 0
                    }

                    DD = DD %>% sort()
                    if (max_norm == TRUE) {
                        DD = DD / max(DD, na.rm = TRUE)
                    }

                    ## max(DD) can be zero
                    if (any(is.na(DD)) | any(is.infinite(DD))) {
                        DD[which(is.na(DD))] = 0
                        DD[which(is.infinite(DD))] = 0
                    }

                    if (size_different) {
                        Position <- seq(from = 0, to = 1, length = length(DD))
                        Points <- tibble(Position = Position, DD = DD)

                        Position_c <- seq(from = 0, to = 1, length = Length)
                        Points_c <- apply(Points, 2, function(u) spline(Position, u, xout = Position_c)$y) %>% as_tibble()

                        # DD_combined = c(DD_combined, Points_c$DD)
                        DD_combined[[DD_kind_name]] = Points_c$DD
                    } else {
                        DD_combined[[DD_kind_name]] = DD
                    }


                }


                # DD_difference_each <- {}
                DD_difference_each <- rep(NA, length(DD_kind))
                for (DD_kind_name in seq_along(DD_kind)) {
                    DD_1 <- DD_target_combined[[DD_kind_name]]
                    DD_2 <- DD_combined[[DD_kind_name]]

                    DD_1_vs_2 <- ((DD_1 - DD_2)^2) #%>% sum() %>% sqrt()
                    
                    if (max(DD_1_vs_2, na.rm = TRUE) != 0) {
                        DD_1_vs_2 = (sum(DD_1_vs_2) / max(DD_1_vs_2, na.rm = TRUE)) * DD_weight[DD_kind_name]
                    } else {
                        DD_1_vs_2 = 0
                    }

                    DD_difference_each[DD_kind_name] = DD_1_vs_2
                }

                DD_difference <- sum(DD_difference_each, na.rm = TRUE) %>% sqrt()

                D_DD[net] = DD_difference

            }

            return_vector <- D_DD










        } else if (net_kind == "list") {

            ## Give weights of one if missing
            if (ncol(target) == 2) {
                target = cbind(target, 1)
            }

            DD_target <- switch(
                DD_kind,
                "out" = {
                            DD_vector <- {}
                            for (n in 1:net_size) {
                                relevant_edges <- target[,1] == n
                                DD_node <- target[relevant_edges, 3] %>% sum()
                                DD_vector = c(DD_vector, DD_node)
                            }
                            DD_vector
                        },
                "in" = {
                            DD_vector <- {}
                            for (n in 1:net_size) {
                                relevant_edges <- target[,2] == n
                                DD_node <- target[relevant_edges, 3] %>% sum()
                                DD_vector = c(DD_vector, DD_node)
                            }
                            DD_vector
                        },
                "undirected" = {
                            DD_vector <- {}
                            for (n in 1:net_size) {
                                relevant_edges <- (target[,1] == n) | (target[,2] == n)
                                intersection_eges <- (target[,1] == n) & (target[,2] == n)
                                DD_node <- sum(target[relevant_edges, 3]) - sum(target[intersection_eges, 3])
                                DD_vector = c(DD_vector, DD_node)
                            }
                            DD_vector
                        },
                "all" = {
                            DD_vector <- {}
                            for (n in 1:net_size) {
                                relevant_edges_out <- target[,1] == n
                                DD_node_out <- target[relevant_edges_out, 3] %>% sum()

                                relevent_edges_in <- target[,2] == n
                                DD_node_int <- target[relevent_edges_in, 3] %>% sum()

                                DD_vector = c(DD_vector, DD_node_out, DD_node_in)
                            }
                            DD_vector
                        },
                stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
            )

            DD_target = DD_target %>% sort()
            if (max_norm == TRUE) {
                DD_target = DD_target / max(DD_target)
            }

            Length <- length(DD_target)


            if (verbose == TRUE) {
                print(paste0("Comparing the degree distributions of ", length(networks), " networks."))
            }


            D_DD <- rep(NA, times = length(networks))
            for (net in seq_along(networks)) {

                comparison_net <- networks[[net]]

                ## Some networks from make_Systematic() will have zero edges
                if (nrow(comparison_net) == 0) {
                    DD <- rep(0, net_size)
                } else {
                    nodes <- c(comparison_net[,1], comparison_net[,2]) %>% unique()

                    ## Give weights of one if missing
                    if (ncol(comparison_net) == 2) {
                        comparison_net = cbind(comparison_net, 1)
                    }

                    ## NOTE: This block assumes the comparison networks are of the same size as the target network
                    DD <- switch(
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
                                        DD_node_int <- comparison_net[relevent_edges_in, 3] %>% sum()

                                        DD_vector = c(DD_vector, DD_node_out, DD_node_in)
                                    }
                                    DD_vector
                                },
                        stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                    )

                    DD = DD %>% sort()
                    if (max_norm == TRUE) {
                        DD = DD / max(DD)
                    }
                }

                # if (length(DD_target) != length(DD)) {
                #     size_different = TRUE
                # } else {
                #     size_different = FALSE
                # }

                if (size_different) {
                    Position <- seq(from = 0, to = 1, length = length(DD))
                    Points <- tibble(Position = Position, DD = DD)

                    Position_c <- seq(from = 0, to = 1, length = Length)
                    Points_c <- apply(Points, 2, function(u) spline(Position, u, xout = Position_c)$y) %>% as_tibble()

                    D_DD[net] = ((DD_target - Points_c$DD)^2) %>% sum() %>% sqrt()
                } else {
                    D_DD[net] = ((DD_target - DD)^2) %>% sum() %>% sqrt()
                }



            }

            return_vector <- D_DD






        } else {
            stop("Unknown network kind. Must be `list` or `matrix`.")
        }

    } else {
        stop("Method not supported.")
    }

    return(return_vector)
}
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

compare <- function(networks, method, net_size, cause_orientation = "row", DD_kind = "all", DD_resize = "smaller", max_norm = FALSE, cores = 1, diffusion_sampling = 2, diffusion_limit = 10, verbose = FALSE) {
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
                        if (verbose == TRUE) { print(c(net_1, net_2)) }

                        alignment <- netcom::align(networks[[net_1]], 
                                                    networks[[net_2]], 
                                                    base = diffusion_sampling,
                                                    max_duration = diffusion_limit,
                                                    characterization = "entropy", 
                                                    normalization = FALSE)$score
                        D_netcom[net_1, net_2] = alignment
                        D_netcom[net_2, net_1] = alignment
                    }
                }
            }
            diag(D_netcom) = 0

        } else { ## Cores > 1
            ## Use parallelization
            cluster <- parallel::makeCluster(cores, outfile = "")
            doParallel::registerDoParallel(cluster)

            ## Network Alignment
            D_netcom <- foreach (net_1 = 1:length(networks), .combine = rbind, .packages = c("tibble", "dplyr", "netcom")) %dopar% {        
                j_output <- rep(NA, length(networks))
                for (net_2 in 1:length(networks)) {
                     if (verbose == TRUE) { print(c(net_1, net_2)) }
                
                    j_output[net_2] <- netcom::align(networks[[net_1]], 
                                                networks[[net_2]], 
                                                base = diffusion_sampling,
                                                max_duration = diffusion_limit, 
                                                characterization = "entropy", 
                                                normalization = FALSE)$score
                }
                
                return(j_output) 
            }

            stopCluster(cluster)
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

            for (net_1 in seq_along(networks)) {
                for (net_2 in seq_along(networks)) {
                    if(net_2 > net_1) {
                        if (verbose == TRUE) { print(c(net_1, net_2)) }


                        if (DD_resize == "smaller") {
                            net_size <- min(nrow(networks[[net_1]]), nrow(networks[[net_2]]))
                        } else if (DD_resize == "larger") {
                            net_size <- max(nrow(networks[[net_1]]), nrow(networks[[net_2]]))
                        } else {
                            stop("DD_resize parameter unknown.")
                        }



                        DD_1_combined <- {}
                        DD_2_combined <- {}

                        for (DD_kind_name in DD_kind) {

                            ## First network's Degree Distribution
                            DD_1 <- switch(
                                DD_kind_name,
                                "out" = networks[[net_1]] %>% rowSums() %>% as.numeric(),
                                "in" = networks[[net_1]] %>% colSums() %>% as.numeric(),
                                "undirected" = as.numeric(rowSums(networks[[net_1]])) + as.numeric(colSums(networks[[net_1]])),
                                "all" = 0.5 * (rbind( rowSums(networks[[net_1]]) + colSums(networks[[net_1]]) )),
                                "entropy_out" = networks[[net_1]] %>% t() %>% vegan::diversity() %>% as.numeric(),
                                "entropy_in" = networks[[net_1]] %>% vegan::diversity() %>% as.numeric(),
                                "entropy_all" = c(networks[[net_1]] %>% t() %>% vegan::diversity() %>% as.numeric(), networks[[net_1]] %>% vegan::diversity() %>% as.numeric()),
                                "clustering_coefficient" = igraph::transitivity(igraph::graph_from_adjacency_matrix(networks[[net_1]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA), type = "weighted"),
                                "cohesion" = igraph::vertex_connectivity(igraph::graph_from_adjacency_matrix(networks[[net_1]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)),
                                "spectral_decomposition" = igraph::embed_adjacency_matrix(igraph::graph_from_adjacency_matrix(networks[[net_1]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA), no = nrow(networks[[net_1]])-1)$X %>% rowSums(),
                                stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                            )

                            ## Some network properties return NA or NaN or Inf on disconnected networks
                            if (any(is.na(DD_1)) | any(is.infinite(DD_1))) {
                                DD_1[which(is.na(DD_1))] = 0
                                DD_1[which(is.infinite(DD_1))] = 0
                            }

                            DD_1 = DD_1 %>% sort()
                            if (max_norm == TRUE) {
                                DD_1 = DD_1 / max(DD_1, na.rm = TRUE)
                            }

                            ## Some network properties return NA or NaN or Inf on disconnected networks
                            if (any(is.na(DD_1)) | any(is.infinite(DD_1))) {
                                DD_1[which(is.na(DD_1))] = 0
                                DD_1[which(is.infinite(DD_1))] = 0
                            }

                            Position_1 <- seq(from = 0, to = 1, length = length(DD_1))
                            Points_1 <- tibble(Position = Position_1, DD = DD_1)

                            Position_c_1 <- seq(from = 0, to = 1, length = net_size)
                            Points_c_1 <- apply(Points_1, 2, function(u) spline(Position_1, u, xout = Position_c_1)$y) %>% as_tibble()

                        ## Second network's Degree Distribution
                            DD_2 <- switch(
                                DD_kind_name,
                                "out" = networks[[net_2]] %>% rowSums() %>% as.numeric(),
                                "in" = networks[[net_2]] %>% colSums() %>% as.numeric(),
                                "undirected" = as.numeric(rowSums(networks[[net_2]])) + as.numeric(colSums(networks[[net_2]])),
                                "all" = 0.5 * (rbind( rowSums(networks[[net_2]]) + colSums(networks[[net_2]]) )),
                                "entropy_out" = networks[[net_2]] %>% t() %>% vegan::diversity() %>% as.numeric(),
                                "entropy_in" = networks[[net_2]] %>% vegan::diversity() %>% as.numeric(),
                                "entropy_all" = c(networks[[net_2]] %>% t() %>% vegan::diversity() %>% as.numeric(), networks[[net_2]] %>% vegan::diversity() %>% as.numeric()),
                                "clustering_coefficient" = igraph::transitivity(igraph::graph_from_adjacency_matrix(networks[[net_2]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA), type = "weighted"),
                                "cohesion" = igraph::vertex_connectivity(igraph::graph_from_adjacency_matrix(networks[[net_2]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)),
                                "spectral_decomposition" = igraph::embed_adjacency_matrix(igraph::graph_from_adjacency_matrix(networks[[net_2]], mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA), no = nrow(networks[[net_2]])-1)$X %>% rowSums(),
                                stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                            )

                            ## Some network properties return NA or NaN or Inf on disconnected networks
                            if (any(is.na(DD_2)) | any(is.infinite(DD_2))) {
                                DD_2[which(is.na(DD_2))] = 0
                                DD_2[which(is.infinite(DD_2))] = 0
                            }
                            
                            DD_2 = DD_2 %>% sort()
                            if (max_norm == TRUE) {
                                DD_2 = DD_2 / max(DD_2, na.rm = TRUE)
                            }

                            ## Some network properties return NA or NaN or Inf on disconnected networks
                            if (any(is.na(DD_2)) | any(is.infinite(DD_2))) {
                                DD_2[which(is.na(DD_2))] = 0
                                DD_2[which(is.infinite(DD_2))] = 0
                            }

                            Position_2 <- seq(from = 0, to = 1, length = length(DD_2))
                            Points_2 <- tibble(Position = Position_2, DD = DD_2)

                            Position_c_2 <- seq(from = 0, to = 1, length = net_size)
                            Points_c_2 <- apply(Points_2, 2, function(u) spline(Position_2, u, xout = Position_c_2)$y) %>% as_tibble()

                            ## Comapre the Degree Distributions of the two networks
                            DD_1 = Points_c_1$DD
                            DD_2 = Points_c_2$DD

                            DD_1_combined = c(DD_1_combined, DD_1)
                            DD_2_combined = c(DD_2_combined, DD_2)
                        }

                        DD_difference <- ((DD_1 - DD_2)^2) %>% sum() %>% sqrt()

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
                        if (verbose == TRUE) { print(c(net_1, net_2)) }

                        ## Skip because edge lists do not inherently contain the number of nodes (zero degree nodes are missing)
                        ## Assume all networks have size net_size
                        # if (DD_resize == "smaller")
                        #     net_size <- min(nrow(networks[[net_1]]), nrow(networks[[net_2]]))
                        # else if (DD_resize == "larger") {
                        #     net_size <- max(nrow(networks[[net_1]]), nrow(networks[[net_2]]))
                        # } else {
                        #     stop("DD_resize parameter unknown.")
                        # }

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
                                                DD_node_int <- comparison_net[relevent_edges_in, 3] %>% sum()

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


                        Position_1 <- seq(from = 0, to = 1, length = length(DD_1))
                        Points_1 <- tibble(Position = Position_1, DD = DD_1)

                        Position_c_1 <- seq(from = 0, to = 1, length = net_size)
                        Points_c_1 <- apply(Points_1, 2, function(u) spline(Position_1, u, xout = Position_c_1)$y) %>% as_tibble()







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
                                                DD_node_int <- comparison_net[relevent_edges_in, 3] %>% sum()

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


                        Position_2 <- seq(from = 0, to = 1, length = length(DD_2))
                        Points_2 <- tibble(Position = Position_2, DD = DD_2)

                        Position_c_2 <- seq(from = 0, to = 1, length = net_size)
                        Points_c_2 <- apply(Points_2, 2, function(u) spline(Position_2, u, xout = Position_c_2)$y) %>% as_tibble()



                        ## Comapre the Degree Distributions of the two networks
                        DD_1 <- Points_c_1$DD
                        DD_2 <- Points_c_2$DD
                        DD_difference <- ((DD_1 - DD_2)^2) %>% sum() %>% sqrt()

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

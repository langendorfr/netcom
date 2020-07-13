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

compare_Target <- function(target, networks, net_size, net_kind, method, cause_orientation = "row", DD_kind = "out", max_norm = FALSE, cores = 1, verbose = FALSE) {
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

            DD_target <- switch(
                DD_kind,
                "out" = target %>% rowSums() %>% as.numeric(),
                "in" = target %>% colSums() %>% as.numeric(),
                "undirected" = as.numeric(rowSums(target)) + as.numeric(colSums(target)),
                stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
            )

            DD_target = DD_target %>% sort()
            if (max_norm == TRUE) {
                DD_target = DD_target / max(DD_target)
            }

            Length <- length(DD_target)

            D_DD <- rep(NA, times = length(networks))
            for (net in seq_along(networks)) {
                if (verbose == TRUE) {
                    print(paste0("Comparing the degree distribution of ", net, " of ", length(networks), " networks."))
                }

                DD <- switch(
                    DD_kind,
                    "out" = networks[[net]] %>% rowSums() %>% as.numeric(),
                    "in" = networks[[net]] %>% colSums() %>% as.numeric(),
                    "undirected" = as.numeric(rowSums(networks[[net]])) + as.numeric(colSums(networks[[net]])),
                    stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                )

                DD = DD %>% sort()
                if (max_norm == TRUE) {
                    DD = DD / max(DD)
                }

                Position <- seq(from = 0, to = 1, length = length(DD))
                Points <- tibble(Position = Position, DD = DD)

                Position_c <- seq(from = 0, to = 1, length = Length)
                Points_c <- apply(Points, 2, function(u) spline(Position, u, xout = Position_c)$y) %>% as_tibble()

                D_DD[net] = ((DD_target - Points_c$DD)^2) %>% sum() %>% sqrt()

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
                stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
            )

            DD_target = DD_target %>% sort()
            if (max_norm == TRUE) {
                DD_target = DD_target / max(DD_target)
            }

            Length <- length(DD_target)





            D_DD <- rep(NA, times = length(networks))
            for (net in seq_along(networks)) {
                if (verbose == TRUE) {
                    print(paste0("Comparing the degree distribution of ", net, " of ", length(networks), " networks."))
                }

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
                        stop("Kind of Degree Distribution not supported. Check the DD_kind parameter.")
                    )

                    DD = DD %>% sort()
                    if (max_norm == TRUE) {
                        DD = DD / max(DD)
                    }
                }


                Position <- seq(from = 0, to = 1, length = length(DD))
                Points <- tibble(Position = Position, DD = DD)

                Position_c <- seq(from = 0, to = 1, length = Length)
                Points_c <- apply(Points, 2, function(u) spline(Position, u, xout = Position_c)$y) %>% as_tibble()

                D_DD[net] = ((DD_target - Points_c$DD)^2) %>% sum() %>% sqrt()

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

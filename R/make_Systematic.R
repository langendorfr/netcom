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

make_Systematic <- function(net_size, net_kind = "list", resolution = 100, resolution_min = 0.01, resolution_max = 0.99, reps = 3, processes = c("ER", "PA", "DD", "SW", "NM"), power_max = 5, cores = 6, directed = TRUE, verbose = TRUE) {
    ## Libraries

    ### Main Body ---
    networks <- list()
    parameters <- tibble(Process = character(),
                        Parameter_Name = character(),
                        Parameter_Value = numeric())
    master_par_systematic <- seq(from = resolution_min, to = resolution_max, length.out = resolution)

    counter <- 0
    for (p in seq_along(processes)) {
        if (verbose == TRUE) {
            print(paste0("Generating ", processes[p], " systems."))
        }
        for (i in 1:resolution) {
            for (r in 1:reps) {
                counter = counter + 1

                if (verbose == TRUE) {
                    # print(c(processes[p], i, r))
                }

                if (processes[p] == "ER") {
                    p_ER <- master_par_systematic[i]
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

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "p_ER", Parameter_Value = p_ER)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "PA") {
                    power_PA <- power_max * master_par_systematic[i]
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

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "power_PA", Parameter_Value = power_PA)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "DD") {
                    divergence_DD <- master_par_systematic[i]
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

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "divergence_DD", Parameter_Value = divergence_DD)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "SW") {
                    rewire_SW <- master_par_systematic[i]
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

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "rewire_SW", Parameter_Value = rewire_SW)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "NM") {
                    connectance_NM <- master_par_systematic[i]
                    niches <- runif(net_size) # %>% sort()
                    net <- make_NM(size = net_size,
                                   net_kind = net_kind,
                                   niches = niches, 
                                   connectance = connectance_NM, 
                                   directed = directed)
                    networks[[counter]] <- net

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "connectance_NM", Parameter_Value = connectance_NM)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)
                    
                } else {
                    stop("An unknown process was included in the simulation.")
                }

            } ## reps within in parameter value
        } ## resolution, through parameter space, within each process
    } ## processes

    # ## Test for symmetry in matrices (only matters if undirected = TRUE)
    # sapply(networks, isSymmetric) %>% all()

    return_list <- list(
        networks = networks,
        parameters = parameters
    )

    return(return_list)
}
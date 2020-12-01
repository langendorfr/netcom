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

make_Systematic <- function(net_size, neighborhood, directed, net_kind = "matrix", resolution = 100, resolution_min = 0.01, resolution_max = 0.99, reps = 3, processes = c("ER", "PA", "DM", "SW", "NM"), power_max = 5, connectance_max = 0.5, divergence_max = 0.5, mutation_max = 0.5, cores = 1, verbose = TRUE) {
    ## Libraries

    ### Main Body ---
    if (!(net_kind %in% c("matrix", "list"))) {
        stop("Unknown net_kind. Must be `list` or `matrix`.")
    }

    # if (!(net_kind == "matrix")) {
    #     stop("Unknown net_kind. Must be `matrix`. `list` will be supported in future versions.")
    # }

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
                    directed = TRUE

                    p_ER <- master_par_systematic[i]
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
                        edgelist <- net %>% igraph::as.directed(mode = "mutual") %>% igraph::as_edgelist(names = TRUE)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "p_ER", Parameter_Value = p_ER)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "PA") {
                    directed = TRUE

                    power_PA <- power_max * master_par_systematic[i]
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

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "power_PA", Parameter_Value = power_PA)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "DD") {
                    directed = FALSE

                    divergence_DD <- divergence_max * master_par_systematic[i]

                    net <- make_DD(size = net_size, 
                                   net_kind = net_kind, 
                                   divergence = divergence_DD, 
                                   directed = directed) #FALSE)
 
                    networks[[counter]] = net

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "divergence_DD", Parameter_Value = divergence_DD)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "DM") {
                    directed = FALSE

                    divergence_DM <- divergence_max * master_par_systematic[i]
                    mutation_DM <- mutation_max * master_par_systematic[i]

                    net <- make_DM(size = net_size, 
                                    net_kind = net_kind, 
                                    divergence = divergence_DM,
                                    mutation = 0, #mutation_DM, #0, #0.01,
                                    directed = directed) #FALSE)

                    networks[[counter]] = net

                    ## Note for this divergence_DM = mutation_DM
                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "mutation_DM", Parameter_Value = mutation_DM)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "SW") {
                    directed = FALSE

                    rewire_SW <- master_par_systematic[i]

                    ## SW neighborhood parameter based on net_size if missing
                    if (missing(neighborhood)) {
                        neighborhood = max(1, round(0.1 * net_size))
                    }

                    net <- make_SW(size = net_size, 
                                   net_kind = net_kind, 
                                   rewire = rewire_SW, 
                                   neighborhood = neighborhood, 
                                   directed = directed) #FALSE)

                    networks[[counter]] <- net

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "rewire_SW", Parameter_Value = rewire_SW)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "NM") {
                    directed = TRUE

                    connectance_NM <- connectance_max * master_par_systematic[i]
                    niches <- runif(net_size) # %>% sort()
                    net <- make_NM(size = net_size,
                                   net_kind = net_kind,
                                   niches = niches, 
                                   connectance = connectance_NM, 
                                   directed = directed, #TRUE
                                   grow = TRUE)
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
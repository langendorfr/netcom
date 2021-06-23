#' @title Systematically Make Networks
#'
#' @description Creates a list of networks that systematically spans mechanisms and their respective parameters.
#'
#' @param net_size Number of nodes in the network.
#' 
#' @param neighborhood The range of nodes that form connected communities. Note: This implementation results in overlap of communities.
#' 
#' @param directed Whether the target network is directed. Defaults to TRUE.
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
#' @param processes Defaults to c("ER", "PA", "DD", "SW", "NM"). Vector of process abbreviations. Currently only the default five are supported. Future versions will accept user-defined network-generating functions and associated parameters. ER = Erdos-Renyi random. PA = Preferential Attachment. DD = Duplication and Divergence. SW = Small World. NM = Niche Model.
#' 
#' @param power_max Defaults to 5. The maximum power of attachment in the Preferential Attachment process (PA).
#' 
#' @param connectance_max Defaults to 0.5. The maximum connectance parameter for the Niche Model.
#' 
#' @param divergence_max Defaults to 0.5. The maximum divergence parameter for the Duplication and Divergence/Mutation mechanisms.
#' 
#' @param mutation_max Defaults to 0.5. The maximum mutation parameter for the Duplication and Mutation mechanism.
#' 
#' @param canonical Defautls to FALSE. If TRUE the mechanisms are directed or undirected in accordance with their canonical forms. This negates the value of `directed`.
#' 
#' @param cores = Defaults to 1. The number of cores to run the classification on. When set to 1 parallelization will be ignored.
#' 
#' @param verbose = Defaults to TRUE. Whether to print all messages.
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

make_Systematic_mixture <- function(net_size, neighborhood, directed = TRUE, net_kind = "matrix", resolution = 100, resolution_min = 0.01, resolution_max = 0.99, reps = 3, processes = c("ER", "PA", "DM", "SW", "NM"), power_max = 5, connectance_max = 0.5, divergence_max = 0.5, mutation_max = 0.5, canonical = FALSE, cores = 1, verbose = TRUE) {
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

                # if (verbose == TRUE) { print(c(processes[p], i, r)) }

                if (processes[p] == "ER") {
                    if (canonical) {directed = TRUE}

                    p_ER <- master_par_systematic[i]
                    mat <- make_Mixture(sequence = rep(paste0("g", processes[p], "_", p_ER), net_size),
                                        # p_ER = p_ER,
                                        directed = directed)

                    if (net_kind == "matrix") {
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
                    if (canonical) {directed = TRUE}

                    power_PA <- power_max * master_par_systematic[i]
                    mat <- make_Mixture(sequence = rep(paste0("g", processes[p], "_", power_PA), net_size),
                                        # power_PA = power_PA,
                                        directed = directed)

                    if (net_kind == "matrix") {
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
                    if (canonical) {directed = FALSE}

                    divergence_DD <- divergence_max * master_par_systematic[i]

                    mat <- make_Mixture(sequence = rep(paste0("g", processes[p], "_", divergence_DD), net_size),
                                        # divergence_DD = divergence_DD,
                                        directed = directed)

                    if (net_kind == "matrix") {
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- net %>% igraph::as.directed(mode = "mutual") %>% igraph::as_edgelist(names = TRUE)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "divergence_DD", Parameter_Value = divergence_DD)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "DM") {
                    if (canonical) {directed = FALSE}

                    divergence_DM <- divergence_max * master_par_systematic[i]
                    mutation_DM <- mutation_max * master_par_systematic[i]

                    mat <- make_Mixture(sequence = rep(paste0("g", processes[p], "_", divergence_DM), net_size),
                                        # divergence_DM = divergence_DM,
                                        # mutation_DM = mutation_DM,
                                        directed = directed)

                    if (net_kind == "matrix") {
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- net %>% igraph::as.directed(mode = "mutual") %>% igraph::as_edgelist(names = TRUE)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                    ## Note for this divergence_DM = mutation_DM
                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "mutation_DM", Parameter_Value = mutation_DM)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "SW") {
                    if (canonical) {directed = FALSE}

                    rewire_SW <- master_par_systematic[i]

                    # ## SW neighborhood parameter based on net_size if missing
                    # if (missing(neighborhood)) {
                    #     neighborhood = max(1, round(0.1 * net_size))
                    # }

                    mat <- make_Mixture(sequence = rep(paste0("g", processes[p], "_", rewire_SW), net_size),
                                        # rewire_SW = rewire_SW,
                                        directed = directed)

                    if (net_kind == "matrix") {
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- net %>% igraph::as.directed(mode = "mutual") %>% igraph::as_edgelist(names = TRUE)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

                    parameters_addition <- tibble(Process = processes[p], Parameter_Name = "rewire_SW", Parameter_Value = rewire_SW)
                    parameters = dplyr::bind_rows(parameters, parameters_addition)

                } else if (processes[p] == "NM") {
                    if (canonical) {directed = TRUE}

                    connectance_NM <- connectance_max * master_par_systematic[i]
                    niches <- runif(net_size) # %>% sort()

                    mat <- make_Mixture(sequence = rep(paste0("g", processes[p], "_", connectance_NM), net_size),
                                        niches = niches,
                                        # connectance_NM = connectance_NM,
                                        directed = directed)

                    if (net_kind == "matrix") {
                        networks[[counter]] <- mat

                    } else if (net_kind == "list") {
                        edgelist <- net %>% igraph::as.directed(mode = "mutual") %>% igraph::as_edgelist(names = TRUE)
                        networks[[counter]] = edgelist

                    } else {
                        stop("Unknown network kind. Must be `list` or `matrix`.")
                    }

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
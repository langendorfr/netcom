grow_Mixture_properties_extinction <- function(sequence, disturbance = 0.1, disturbance_kind = "proportion", extinction_step = 5, extinction_probability = 0.5, extinction_kind = "edge", p_ER = 0.5, power_PA = 2, divergence_DD = 0.8, link_DD = 0.2, divergence_DM = 0.8, mutation_DM = 0.4, link_DM = 0.2, rewire_SW = 0.1, connectance_NM = 1, retcon = FALSE, directed = TRUE, force_connected = FALSE) {
    size <- length(sequence)
    matrix <- matrix(0, size, size)

    ## Start with the first two nodes connected
    matrix[1,2] = 1
    matrix[2,1] = 1

    ## Assign niches for NM (the Niche Model) so they will be constant across calls to grow_CM
    niches <- runif(size) %>% sort()

    ## Assign clusters for SW (the Small World model) so they will be constant across calls to grow_SW
    # k <- 3
    # num_clusters <- ceiling(size/k)
    # clusters <- rep(1:num_clusters, each = k)
    # clusters = clusters[1:size]

    ## Structural properties to output
    Ascendency <- rep(NA, size-2)
    ACratio <- rep(NA, size-2)

    x <- 3
    step <- 1
    while (x <= size) {
        # print(x)
        if (sequence[x] == "ER") {
            matrix = grow_ER(matrix, x, p = p_ER, retcon = retcon, directed = directed)
        } else if (sequence[x] == "PA") {
            matrix = grow_PA(matrix, x, power = power_PA, retcon = retcon, directed = directed)
        } else if (sequence[x] == "DD") {
            matrix = grow_DD(matrix, x, divergence = divergence_DD, link = link_DD, directed = directed)
        } else if (sequence[x] == "DM") {
            matrix = grow_DM(matrix, x, divergence = divergence_DM, mutation = mutation_DM, link = link_DM, directed = directed)
        } else if (sequence[x] == "SW") {
            matrix = grow_SW(matrix, x, rewire = rewire_SW, retcon = retcon, directed = directed)
        } else if (sequence[x] == "NM") {
            matrix = grow_NM(matrix, x, connectance = connectance_NM, niches = niches, retcon = retcon, directed = directed) 
        } else {
            print("ERROR: Model not specified.")
        }

        if (disturbance_kind == "number") {
            for (d in 1:disturbance) {
                matrix[1:x, 1:x] = stir_Mixture_nodeLevel(matrix = matrix[1:x, 1:x], 
                                                          sequence = rep(process, x), 
                                                          stirs = disturbance, 
                                                          p_ER = p_ER, 
                                                          power_PA = power_PA, 
                                                          divergence_DM = divergence_DM,
                                                          mutation_DM = mutation_DM,
                                                          link_DM = link_DM, 
                                                          connectance_NM = connectance_NM, 
                                                          rewire_SW = rewire_SW, 
                                                          force_connected = force_connected, 
                                                          directed = directed)
            }
        } else if (disturbance_kind == "proportion") {
            stirs = round(disturbance * x, digits = 0)
            if (stirs > 0) {

                for (s in 1:stirs) {
                    matrix[1:x, 1:x] = stir_Mixture_nodeLevel(matrix = matrix[1:x, 1:x], 
                                                              sequence = rep(process, x), 
                                                              stirs = stirs, 
                                                              p_ER = p_ER, 
                                                              power_PA = power_PA, 
                                                              divergence_DM = divergence_DM,
                                                              mutation_DM = mutation_DM,
                                                              link_DM = link_DM, 
                                                              connectance_NM = connectance_NM, 
                                                              rewire_SW = rewire_SW, 
                                                              force_connected = force_connected, 
                                                              directed = directed)
                }

            }
        } else {
            stop("Unrecognized disturbance_kind. Must be `number` or `proportion`.")
        }


        NI_AI <- NetIndices::AscInd(Flow = matrix[1:x, 1:x])
        Ascendency[x-2] = NI_AI["Total", "Ascendency"]
        ACratio[x-2] = NI_AI["Total", "ACratio"]

        if (step == extinction_step) {
            if (extinction_kind == "node") {
                stop("`node` extinction_kind not yet supported.")
            } else if (extinction_kind == "edge") {
                extinction_edge_ids <- which(matrix[1:x, 1:x] != 0)

                for (edge in seq_along(extinction_edge_ids)) {
                    if (runif(1) <= extinction_probability) {
                        matrix[1:x, 1:x][extinction_edge_ids[edge]] = 0
                    }
                }
            } else {
                stop("Unknown extinction_kind. Must be `node` or `edge`.")
            }
        }

        x = x + 1
        step = step + 1
    }

    return(list(Network = matrix,
                Ascendency = Ascendency, 
                ACratio = ACratio))
}



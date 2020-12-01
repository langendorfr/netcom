grow_Mixture <- function(sequence, disturbance = 0, p_ER = 0.5, power_PA = 2, divergence_DD = 0.8, link_DD = 0.2, divergence_DM = 0.8, mutation_DM = 0.4, link_DM = 0.2, rewire_SW = 0.1, connectance_NM = 1, retcon = FALSE, directed = TRUE, force_connected = FALSE) {
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

    for (x in 3:size) {
        # print(x)
        if (sequence[x] == "ER") {
            matrix = grow_ER(matrix, x, p = p_ER, retcon = retcon, directed = TRUE) #directed)
        } else if (sequence[x] == "PA") {
            matrix = grow_PA(matrix, x, power = power_PA, retcon = retcon, directed = TRUE) #directed)
        } else if (sequence[x] == "DD") {
            matrix = grow_DD(matrix, x, divergence = divergence_DD, link = link_DD, directed = FALSE) #directed)
        } else if (sequence[x] == "DM") {
            matrix = grow_DM(matrix, x, divergence = divergence_DM, mutation = mutation_DM, link = link_DM, directed = FALSE) #directed)
        } else if (sequence[x] == "SW") {
            matrix = grow_SW(matrix, x, rewire = rewire_SW, retcon = retcon, directed = FALSE) #directed)
        } else if (sequence[x] == "NM") {
            matrix = grow_NM(matrix, x, connectance = connectance_NM, niches = niches, retcon = retcon, directed = TRUE) #directed) 
        } else {
            print("ERROR: Model not specified.")
        }

        if (disturbance > 0) {
            for (d in 1:disturbance) {
                stir_Mixture(matrix = matrix, 
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
        }
    }

    return(matrix)
}



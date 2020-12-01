stir_Mixture <- function(matrix, sequence, stirs = 100, p_ER = 0.5, power_PA = 2, divergence_DM = 0.8, mutation_DM = 0.1, link_DM = 0.2, connectance_NM = 0.2, rewire_SW = 0.1, force_connected = FALSE, directed = TRUE) {
    # ## Primary Directory
    # pd <- "/Users/ryan/Windows/Documents/Post UCB/Research/Relativism"
    # setwd(pd)

    # ## Libraries
    # library("tidyverse")
    # library("vegan")
    # library("igraph")

    # ## Custom Functions
    # source("stir_ER.R")
    # source("stir_PA.R")
    # source("stir_DD.R")
    # source("stir_SW.R")
    # # source("grow_CC.R")

    size <- length(sequence)

    sequence_id = seq_along(sequence)

    niches <- runif(size) %>% sort()

    sequence_stir = {}
    for (i in 1:stirs) {
        sequence_stir = c(sequence_stir, sample(x = sequence_id, size = length(sequence_id), replace = FALSE))
    }

    for (s in seq_along(sequence_stir)) {
        x = sequence_stir[s]

        # print(matrix)
        # print(sequence[x])

        if (sequence[x] == "ER") {
            matrix = stir_ER(matrix = matrix, x = x, p = p_ER)
        } else if (sequence[x] == "PA") {
            matrix = stir_PA(matrix = matrix, x = x, power = power_PA)
        } else if (sequence[x] == "DM") {
            matrix = stir_DM(matrix = matrix, x = x, divergence = divergence_DM, mutation = mutation_DM, link = link_DM, force_connected = force_connected)
        } else if (sequence[x] == "NM") {
            matrix = stir_NM(matrix = matrix, x = x, niches = niches, connectance = connectance_NM, directed = directed)
        } else if (sequence[x] == "SW") {
            matrix = stir_SW(matrix = matrix, x = x, rewire = rewire_SW)
        } else {
            print("ERROR: Model not specified.")
        }
    }

    return(matrix)
}

stir_DD <- function(matrix, x, divergence, link, force_connected = FALSE) {
    ids <- (1:ncol(matrix))[-x]

    DD <- function() {
        duplication <- sample(ids, 1)
        matrix[, x] <- matrix[, duplication]
        matrix[x, ] <- matrix[duplication, ]

        ## Connect new node to copied node
        if (runif(1) < link) {
        matrix[x,duplication] <- 1
        matrix[duplication,x] <- 1
        }

        for (i in ids) {
            ##!! Do this differently than grow_DD, which uses `if (matrix[i, x] == 1)`
            if (matrix[x, i] == 1) {
                if (runif(1) < divergence) {
                # matrix[i, x] <- 0
                matrix[x, i] <- 0
                }
            }
        }
        
        return(matrix) ## DD
    }

    matrix = DD()

    if (force_connected) {
        ## Rerun until no entirely disconnected nodes
        ##!! Do this differently than grow_DD, which uses `while (sum(matrix[, x]) == 0 && sum(matrix[x,]) == 0) {`
        while (sum(matrix[x,]) == 0) {
        matrix = DD()
        }
    }

    return(matrix) ## grow_DD
}
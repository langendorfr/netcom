grow_DM <- function(matrix, x, divergence, mutation = 0, link = 0, connected = FALSE, retcon = FALSE, directed = TRUE) {
    w <- x - 1
    DD <- function() {
        duplication <- sample(1:w, 1)
        matrix[, x] <- matrix[, duplication]
        matrix[x, ] <- matrix[duplication, ]

        ## Connect new node to copied node
        if (runif(1) < link) {
        matrix[x, duplication] <- 1
        matrix[duplication, x] <- 1
        }

        for (i in 1:w) {

            if (matrix[x, i] == 1) {
                if (runif(1) <= divergence) {
                    matrix[x, i] <- 0
                }
            } else if (matrix[x, i] == 0) {
                if (runif(1) <= mutation) {
                    matrix[x, i] <- 0
                }
            } else {
                stop("Weighted edge detected. Only binary networks are supported in this release.")
            }

                                        # if (runif(1) <= divergence) {

                                        #     if (matrix[x, i] == 0) {
                                        #         matrix[x, i] <- 1
                                        #     } else if (matrix[x, i] == 1) {
                                        #         matrix[x, i] <- 0
                                        #     }

                                        # }
        }

        if (retcon == TRUE) {
            stop("Retcon functionality missing.")
        }

        if (directed == FALSE) {
            matrix[1:x, x] = matrix[x, 1:x]
        }

        return(matrix) ## DD
    }

    matrix = DD()

    ## Rerun until no entirely disconnected nodes
    if (connected == TRUE) {
        while (sum(matrix[, x]) == 0 && sum(matrix[x,]) == 0) {
            matrix = DD()
        }
    }

    return(matrix) ## grow_DD
}
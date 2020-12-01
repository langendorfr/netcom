grow_SW <- function(matrix, x, rewire, connected = FALSE, retcon = FALSE, directed = TRUE) {
    w <- x-1

    out_degree <- rowSums(matrix[1:w, 1:w])
    in_degree <- colSums(matrix[1:w, 1:w])    

    out_avg <- mean(out_degree) %>% round(digits = 0)
    in_avg <- mean(in_degree) %>% round(digits = 0)

    out_new <- 1 * ((out_avg - out_degree) >= 1)
    in_new <- 1 * ((in_avg - in_degree) >= 1)

    ## Do not allow no links
    if (connected == TRUE) {
        if (max(out_new) == 0) {
            out_id <- sample(seq_along(out_new), size = 1)
            in_id <- sample(seq_along(in_new), size = 1)

            out_new[out_id] = 1
            in_new[in_id] = 1
        }
    }

    ## Add zero because no self loops; diag(matrix) = 0
    matrix[x, 1:x] = c(out_new, 0)
    
    if (retcon == TRUE) {
        matrix[1:x, x] = c(in_new, 0)
    } else {
        matrix[1:x, x] = 0
    }

    for (col in 1:w) {
    
        ## Do not attempt rewiring when all edges exist
        out_edges <- which(matrix[x, 1:w] != 0)
        if (length(out_edges) != w) {
            for (e in out_edges) {
                if (runif(1) <= rewire) {
                    e_possible <- which(matrix[x, 1:w] == 0)
                    e_new <- sample(e_possible, 1)
                    
                    matrix[x, e] = 0
                    matrix[x, e_new] = 1
                }
            }
        }
    }

    if (directed == FALSE) {
        matrix[1:x, x] = matrix[x, 1:x]
    }

    return(matrix)
}
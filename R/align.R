#' @title Dynamic Network Alignment
#'
#' @description Network alignment by comparing the entropies of diffusion kernels simulated on two networks.
#' Functions to take two networks, either matrices or linked lists, and return a node-level alignment between them.
#'
#' @param matrix_1_input The first network being aligned, either as a matrix or linked list. If the two
#'     networks are of different sizes, it will be easier to interpret the output if this is the smaller one.
#'
#' @param matrix_2_input The second network. Should be the same type (matrix or linked list) as matrix_1_input.
#'
#' @param base The base in the series of time steps to sample the diffusion kernels at. If base = 1 every time step
#'     is sampled. If base = 2, only time steps that are powers of 2 are sampled, etc.
#'
#' @param characterization Determines how the diffusion kernels are characterized. Either "entropy" or "Gini".
#'
#' @param normalization Determines if self-loops should be augmented such that edge weights in network.1 and network.2 are
#'     proportional to those in matrix.1.input and matrix.2.input. FALSE by default because this is innapropriate for
#'     unweighted binary/logical networks where edges indicate only the presense of an interaction.
#'
#' @return A list containing 4 pieces:
#'
#' Score = mean of all alignment scores between nodes in both original networks matrix_1_input and matrix_2_input.
#'
#' Alignment = data frame of the nodes in both networks, sorted numerically by the first network (why it helps to make the smaller network the first one), and the corresponding alignment score
#'
#' Score_with_Padding = same as Score but includes the padding nodes in the smaller network, which can be thought of as a size gap penalty for aligning differently sized networks
#'
#' Alignment_with_Padding = same as Alignment but includes the padding nodes in the smaller network
#'
#' @examples
#' NetOne <- matrix(runif(25,0,1), nrow=5, ncol=5)
#' NetTwo <- matrix(runif(25,0,1), nrow=5, ncol=5)
#' NetCom(NetOne, NetTwo)
#' NetCom(NetOne, NetTwo, input = "matrix", base = 2, characterization = "entropy", normalization = FALSE)
#'

NetCom <- function(matrix_1_input, matrix_2_input, input = "matrix", base = 2, characterization = "entropy", normalization = FALSE)
{
  # Ensure inputs are matrices, and if they are linked lists convert them to matrices (NOTE: this assumes the same data type for the two input networks)
  if (input == "list" | dim(matrix_1_input)[1] != dim(matrix_1_input)[2] | dim(matrix_2_input)[1] != dim(matrix_2_input)[2]) {

    # R starts counting at one, not zero
    if (min(matrix_1_input) == 0) {
      matrix_1_input <- as_adjacency_matrix(graph_from_edgelist(as.matrix(matrix_1_input + 1), directed = TRUE), sparse = FALSE)
    } else {
      matrix_1_input <- as_adjacency_matrix(graph_from_edgelist(as.matrix(matrix_1_input), directed = TRUE), sparse = FALSE)
    }

    if (min(matrix_2_input) == 0) {
      matrix_2_input <- as_adjacency_matrix(graph_from_edgelist(as.matrix(matrix_2_input + 1), directed = TRUE), sparse = FALSE)
    } else {
      matrix_2_input <- as_adjacency_matrix(graph_from_edgelist(as.matrix(matrix_2_input), directed = TRUE), sparse = FALSE)
    }

  }

  # Calculate the number of nodes in each network, which are used repeatedly
  matrix_sizes <- c(nrow(matrix_1_input), nrow(matrix_2_input))

  # Preserve relative edge weights (turned off by default, for unweighted binary networks)
  if (normalization == TRUE) {
    matrix_1 <- matrix_1_input
    row_max_1 <- max(Matrix::rowSums(matrix_1_input))
    for (i in 1:matrix_sizes[1]) {
      matrix_1[i, i] <-  matrix_1[i, i] + (row_max_1 - sum(matrix_1_input[i, ]))
    }

    matrix_2 <- matrix_2_input
    row_max_2 <- max(Matrix::rowSums(matrix_2_input))
    for (i in 1:matrix_sizes[2]) {
      matrix_2[i, i] <-  matrix_2[i, i] + (row_max_2 - sum(matrix_2_input[i, ]))
    }
  } else {
    matrix_1 <- matrix_1_input
    matrix_2 <- matrix_2_input
  }

  # Add padding (disconnected nodes)
  if (which.min(matrix_sizes) == 1) {
    matrix_1_padded <- diag(matrix_sizes[2])
    matrix_1_padded[1:matrix_sizes[1], 1:matrix_sizes[1]] <- matrix_1

    # Only the smaller network has padding added, but for readability both get the new name
    matrix_2_padded <- matrix_2_input
  } else {
    matrix_2_padded <- diag(matrix_sizes[1])
    matrix_2_padded[1:matrix_sizes[2], 1:matrix_sizes[2]] <- matrix_2

    # Only the smaller network has padding added, but for readability both get the new name
    matrix_1_padded <- matrix_1_input
  }

  # Normalize matrices into row-stochastic markov processes
  network_1 <- sweep(matrix_1_padded, 1, Matrix::rowSums(matrix_1_padded), FUN = "/")
  network_2 <- sweep(matrix_2_padded, 1, Matrix::rowSums(matrix_2_padded), FUN = "/")

  # Remove NA from rows of sink nodes
  network_1[is.nan(network_1)] <- 0
  network_2[is.nan(network_2)] <- 0

  # Add 1 to the diagonal in empty rows to prevent diffusion out of the system
  network_1_sinks <- which(rowSums(network_1) == 0)
  network_2_sinks <- which(rowSums(network_2) == 0)

  if (length(network_1_sinks) > 0) {
    for (i in 1:length(network_1_sinks)) {
      network_1[network_1_sinks[i], network_1_sinks[i]] <- 1
    }
  }

  if (length(network_2_sinks) > 0) {
    for (j in 1:length(network_2_sinks)) {
      network_2[network_2_sinks[j], network_2_sinks[j]] <- 1
    }
  }

  # Duration of the simulation (note: padded nodes are disconnected and therefore do not change a network's diameter)
  network_1_duration <- length(igraph::get.diameter(igraph::graph.adjacency(network_1, weighted = TRUE, mode = 'directed', diag = TRUE))) - 1   # Subtract 1 because the number of steps in the diameter is one less than the number of nodes in it
  network_2_duration <- length(igraph::get.diameter(igraph::graph.adjacency(network_2, weighted = TRUE, mode = 'directed', diag = TRUE))) - 1   # Subtract 1 because the number of steps in the diameter is one less than the number of nodes in it
  duration_max <- max(2 * network_1_duration, 2 * network_2_duration, 2)                                                                    # The 2 at the end ensures at least two time steps, and the other 2's are for twice the diameter

  # base = 1 is the most conservative (sample every time step) but also the most expensive computationally
  if (missing(base) | base == 1) {
    base <- 1
    kernel_sampling <- 1:duration_max
  } else {
    kernel_sampling <- base^(0:ceiling(log(duration_max, base)))
  }

  for (network in 1:2) {

    if (network == 1) {
      # Matrix of the diffusion kernel statistic (e.g. normalized entropy, Gini coefficient) through time
      network_1_output <- matrix(NA, nrow = max(matrix_sizes), ncol = length(kernel_sampling))

      # Characterize the first time step
      network_1_diffusion <- network_1 # Dummy version of network_1 that will hold the diffusions

      if (characterization == "entropy") {
        network_1_output[, 1] <- c(vegan::diversity(network_1_diffusion[1:matrix_sizes[1], 1:matrix_sizes[1]]) / vegan::diversity(rep(1, matrix_sizes[1])), rep(0, max(matrix_sizes) - matrix_sizes[1]))
      }
      if (characterization == "Gini") {
        network_1_output[, 1] <- c(Gini(network_1_diffusion[1:matrix_sizes[1], 1:matrix_sizes[1]]), rep(0, max(matrix_sizes) - matrix_sizes[1]))
      }

      # Characterize the remaining time steps
      for (t in 2:length(kernel_sampling)) {
        network_1_diffusion <- network_1_diffusion %*% (network_1 %^% (kernel_sampling[t] - kernel_sampling[t - 1])) # Speeds up the algorithm by eliminating the need to recalculate previous time steps

        # For those interested in speeding things up, the following commented-out line proved to be slower
        # network_1_diffusion <- abs(Re(eig$vectors %*% diag(eig$values ^ t) %*% inverse))

        if (characterization == "entropy") {
          network_1_output[, t] <- c(vegan::diversity(network_1_diffusion[1:matrix_sizes[1], 1:matrix_sizes[1]]) / vegan::diversity(rep(1, matrix_sizes[1])), rep(0, max(matrix_sizes) - matrix_sizes[1]))
        }
        if (characterization == "Gini") {
          network_1_output[, t] <- c(Gini(network_1_diffusion[1:matrix_sizes[1], 1:matrix_sizes[1]]), rep(0, max(matrix_sizes) - matrix_sizes[1]))
        }
      }

    } else {
      # Identical process, but for the second network
      network_2_output <- matrix(NA, nrow = max(matrix_sizes), ncol = length(kernel_sampling))

      network_2_diffusion <- network_2
      if (characterization == "entropy") {
        network_2_output[, 1] <- c(vegan::diversity(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]]) / vegan::diversity(rep(1, matrix_sizes[2])), rep(0, max(matrix_sizes) - matrix_sizes[2]))
      }
      if (characterization == "Gini") {
        network_2_output[, 1] <- c(Gini(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]]), rep(0, max(matrix_sizes) - matrix_sizes[2]))
      }

      for (t in 2:length(kernel_sampling)) {
        network_2_diffusion <- network_2_diffusion %*% (network_2 %^% (kernel_sampling[t] - kernel_sampling[t - 1]))
        if (characterization == "entropy") {
          network_2_output[, t] <- c(vegan::diversity(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]]) / vegan::diversity(rep(1, matrix_sizes[2])), rep(0, max(matrix_sizes) - matrix_sizes[2]))
        }
        if (characterization == "Gini") {
          network_2_output[, t] <- c(Gini(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]]), rep(0, max(matrix_sizes) - matrix_sizes[2]))
        }

      }
    }

  }

  # Calculate all pairwise distances between the rows of the characterized diffusion over time matrices (network.1.output and network.2.output) normalized by the number of time steps sampled
  cost_matrix <- as.matrix(suppressWarnings(pdist::pdist(network_1_output, network_2_output))) / length(kernel_sampling) # Warnings are suppressed to prevent being notified when network.1.output = network.2.output

  # Find the optimal assignment using the Hungarian algorithm
  assignment <- clue::solve_LSAP(cost_matrix, maximum = FALSE)
  mapping <- cbind(seq_along(assignment), assignment)

  # The 'padding' versions include the nodes added to the smaller network (if matrix.sizes[1] != matrix.sizes[2]) as a size penalty
  alignment_padding <- cbind(mapping, cost_matrix[mapping])
  colnames(alignment_padding) <- c("Network_1", "Network_2", "Score")

  alignment <- alignment_padding[alignment_padding[, 1] <= matrix_sizes[1] & alignment_padding[, 2] <= matrix_sizes[2], ]
  colnames(alignment) <- c("Network_1", "Network_2", "Score")

  alignment_score_padding <- mean(alignment_padding[, 3]) # Taking the mean standardizes the alignment score across alignments between networks of varying sizes
  alignment_score <- mean(alignment[, 3]) # Ignore the padding nodes, which act as a penalty for aligning networks of different sizes

  output <- list("Score" = alignment_score, "Alignment" = alignment, "Score_with_Padding" = alignment_score_padding, "Alignment_with_Padding" = alignment_padding)


  return(output)
}

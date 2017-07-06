#' @title Dynamic Network Alignment
#'
#' @description Network alignment by comparing the entropies of diffusion kernels simulated on two networks.
#' \code{align} takes two networks, either matrices or linked lists, and returns a node-level alignment between them.
#'
#' @param network_1_input The first network being aligned, either as a matrix or linked list. If the two
#'     networks are of different sizes, it will be easier to interpret the output if this is the smaller one.
#'
#' @param network_2_input The second network. Should be the same type (matrix or linked list) as network_1_input.
#'
#' @param input Defaults to "matrix". Can be set to "list" if the two networks are stored as linked lists.
#'
#' @param base Defaults to 1. The base in the series of time steps to sample the diffusion kernels at. If base = 1 every time step
#'     is sampled. If base = 2, only time steps that are powers of 2 are sampled, etc. Larger values place more emphasis on 
#'     earlier time steps. This can be helpful if the diffusion kernel quickly converges to an equilibrium, and also
#'     runs faster. 
#'     
#' @param max_duration Defaults to twice the diameter of the larger network. Sets the number of time steps to allow the diffusion
#'     kernel to spread for, which is the smallest power of base that is at least as large as max_duration.     
#'
#' @param characterization Defaults to "entropy". Determines how the diffusion kernels are characterized. Either "entropy" or "gini".
#'
#' @param normalization Defaults to FALSE. Determines if self-loops should be augmented such that edge weights are
#'     proportional to those in network_1_input and network_2_input. FALSE by default because this is innapropriate for
#'     unweighted binary/logical networks where edges indicate only the presense of an interaction.
#'     
#' @param unit_test Defaults to FALSE. Saves the following intermediate steps to help with general troubleshooting: post-processing matrix
#'     representations of both networks, time steps at which the diffusion kernels were sampled, the diffusion kernels at those time steps,
#'     the characterizations of the diffusion kernels at those time steps, and the cost matrix fed into the Hungarian algorithm where the 
#'     ij element is the difference between the characterization-over-time curves for node i in the first network and node j in the second network.
#'
#' @details Network alignment pairs nodes between two networks so as to maximize similarities in their edge structures. 
#'     This allows information from well-studied systems to be used in poorly studied ones, such as to identify
#'     unknown protein functions or ecosystems that will respond similarly to a given disturbance. Most network alignment
#'     algorithms focus on measures of topological overlap between edges of the two networks. The method implemented here 
#'     compares nodes using the predictability of dynamics originating from each node in each network. Consider network alignment 
#'     as trying to compare two hypothetical cities of houses connected by roads. The approach implemented here is to pairwise 
#'     compare each house with those in the other city by creating a house-specific signature. This is accomplished by quantifying 
#'     the predictability of the location of a person at various times after they left their house, assuming they were moving randomly. 
#'     This predictability across all houses captures much of the way each city is organized and functions. \code{align} 
#'     uses this conceptual rationale to align two networks, with nodes as houses, edges as roads, and random diffusion representing people leaving 
#'     their houses and walking around the city to other houses. The mechanics of this, which are conceptually akin to flow 
#'     algorithms and Laplacian dynamics, can be analytically expressed as a Markov chain raised to successive powers which are 
#'     the durations of diffusion.
#'
#' @return 
#' \item{score}{Mean of all alignment scores between nodes in both original networks network_1_input and network_2_input.}
#' \item{alignment}{Data frame of the nodes in both networks, sorted numerically by the first network (why it helps to make the smaller network the first one), and the corresponding alignment score.}
#' \item{score_with_padding}{Same as score but includes the padding nodes in the smaller network, which can be thought of as a size gap penalty for aligning differently sized networks. Only included if the input networks are different sizes.}
#' \item{alignment_with_padding}{Same as alignment but includes the padding nodes in the smaller network. Only included if the input networks are different sizes.}
#' 
#' @examples
#' net_one <- matrix(runif(25,0,1), nrow=5, ncol=5)
#' net_two <- matrix(runif(25,0,1), nrow=5, ncol=5)
#' align(net_one, net_two)
#' align(net_one, net_two, base = 1, characterization = "gini", normalization = TRUE)
#' 
#' @export

align <- function(network_1_input, network_2_input, input = "matrix", base = 2, max_duration, characterization = "entropy", normalization = FALSE, unit_test = FALSE)
{
  # Prevents the NOTE "no visible binding for global variable" which arises because these variables are set in a loop with their names deriving from the loop index 
  network_1_diffusion_2 = network_1_diffusion_4 = network_1_diffusion_8 = NULL
  network_2_diffusion_2 = network_2_diffusion_4 = network_2_diffusion_8 = NULL
  network_1_output_2 = network_1_output_4 = network_1_output_8 = NULL
  network_2_output_2 = network_2_output_4 = network_2_output_8 = NULL
  
  # Check if inputs are square matrices. If not, they are linked lists which need to be converted to 
  # their respective matrix representations. (NOTE: this assumes the same data type for the two input networks)
  if (input == "list" | dim(network_1_input)[1] != dim(network_1_input)[2] | dim(network_2_input)[1] != dim(network_2_input)[2]) {

    # R starts counting at one, not zero
    if (min(network_1_input) == 0) {
      matrix_1_input <- igraph::as_adjacency_matrix(igraph::graph_from_edgelist(as.matrix(network_1_input + 1), directed = TRUE), sparse = FALSE)
    } else {
      matrix_1_input <- igraph::as_adjacency_matrix(igraph::graph_from_edgelist(as.matrix(network_1_input), directed = TRUE), sparse = FALSE)
    }

    if (min(network_2_input) == 0) {
      matrix_2_input <- igraph::as_adjacency_matrix(igraph::graph_from_edgelist(as.matrix(network_2_input + 1), directed = TRUE), sparse = FALSE)
    } else {
      matrix_2_input <- igraph::as_adjacency_matrix(igraph::graph_from_edgelist(as.matrix(network_2_input), directed = TRUE), sparse = FALSE)
    }

  } else {
    
    # When the input networks are matrices
    matrix_1_input <- network_1_input
    matrix_2_input <- network_2_input
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
  if (missing(max_duration)) {
    network_1_duration <- length(igraph::get.diameter(igraph::graph.adjacency(network_1, weighted = TRUE, mode = 'directed', diag = TRUE))) - 1   # Subtract 1 because the number of steps in the diameter is one less than the number of nodes in it
    network_2_duration <- length(igraph::get.diameter(igraph::graph.adjacency(network_2, weighted = TRUE, mode = 'directed', diag = TRUE))) - 1   # Subtract 1 because the number of steps in the diameter is one less than the number of nodes in it
    duration_upper_bound <- max(2 * network_1_duration, 2 * network_2_duration, 2)                                                                    # The 2 at the end ensures at least two time steps, and the other 2's are for twice the diameter
  } else {
    duration_upper_bound <- max_duration
  }
  
  # base = 1 is the most conservative (sample every time step) but also the most expensive computationally
  if (missing(base) | base == 1) {
    base <- 1
    kernel_sampling <- 1:duration_upper_bound
  } else {
    kernel_sampling <- base^(0:ceiling(log(duration_upper_bound, base)))
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
      if (characterization == "gini") {
        network_1_output[, 1] <- c(gini(network_1_diffusion[1:matrix_sizes[1], 1:matrix_sizes[1]], byrow = TRUE), rep(0, max(matrix_sizes) - matrix_sizes[1]))
      }
      
      # unit_test functionality
      if (unit_test == TRUE) {
        network_1_diffusion_1 <- network_1_diffusion
        network_1_output_1 <- network_1_output[, 1]
      }

      # Characterize the remaining time steps
      for (t in 2:length(kernel_sampling)) {
        network_1_diffusion <- expm::"%^%"(network_1, kernel_sampling[t]) 

        # For those interested in speeding things up, the following commented-out line proved to be slower
        # network_1_diffusion <- abs(Re(eig$vectors %*% diag(eig$values ^ t) %*% inverse))
        
        if (characterization == "entropy") {
          network_1_output[, t] <- c(vegan::diversity(network_1_diffusion[1:matrix_sizes[1], 1:matrix_sizes[1]]) / vegan::diversity(rep(1, matrix_sizes[1])), rep(0, max(matrix_sizes) - matrix_sizes[1]))
        }
        if (characterization == "gini") {
          network_1_output[, t] <- c(gini(network_1_diffusion[1:matrix_sizes[1], 1:matrix_sizes[1]], byrow = TRUE), rep(0, max(matrix_sizes) - matrix_sizes[1]))
        }
        
        # unit_test functionality
        if (unit_test == TRUE) {
          assign(paste("network_1_diffusion", kernel_sampling[t], sep = "_"), network_1_diffusion)
          assign(paste("network_1_output", kernel_sampling[t], sep = "_"), network_1_output[, t])          
        }
      }

    } else {
      # Identical process, but for the second network
      network_2_output <- matrix(NA, nrow = max(matrix_sizes), ncol = length(kernel_sampling))

      network_2_diffusion <- network_2
      
      if (characterization == "entropy") {
        network_2_output[, 1] <- c(vegan::diversity(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]]) / vegan::diversity(rep(1, matrix_sizes[2])), rep(0, max(matrix_sizes) - matrix_sizes[2]))
      }
      if (characterization == "gini") {
        network_2_output[, 1] <- c(gini(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]], byrow = TRUE), rep(0, max(matrix_sizes) - matrix_sizes[2]))
      }
      
      # unit_test functionality
      if (unit_test == TRUE) {
        network_2_diffusion_1 <- network_2_diffusion
        network_2_output_1 <- network_2_output[, 1]
      }

      for (t in 2:length(kernel_sampling)) {
        network_2_diffusion <- expm::"%^%"(network_2, kernel_sampling[t])
        
        if (characterization == "entropy") {
          network_2_output[, t] <- c(vegan::diversity(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]]) / vegan::diversity(rep(1, matrix_sizes[2])), rep(0, max(matrix_sizes) - matrix_sizes[2]))
        }
        if (characterization == "gini") {
          network_2_output[, t] <- c(gini(network_2_diffusion[1:matrix_sizes[2], 1:matrix_sizes[2]], byrow = TRUE), rep(0, max(matrix_sizes) - matrix_sizes[2]))
        }
        
        # unit_test functionality
        if (unit_test == TRUE) {
          assign(paste("network_2_diffusion", kernel_sampling[t], sep = "_"), network_2_diffusion)
          assign(paste("network_2_output", kernel_sampling[t], sep = "_"), network_2_output[, t])
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

  if (unit_test == TRUE) {
    output <- list("score" = alignment_score, 
                   "alignment" = alignment, 
                   "score_with_padding" = alignment_score_padding, 
                   "alignment_with_padding" = alignment_padding,
                   "network_1" = network_1,
                   "network_2" = network_2,
                   "matrix_1_padded" = matrix_1_padded,
                   "matrix_2_padded" = matrix_2_padded,
                   "kernel_sampling" = kernel_sampling,
                   "network_1_diffusion_1" = network_1_diffusion_1,
                   "network_1_diffusion_2" = network_1_diffusion_2,
                   "network_1_diffusion_4" = network_1_diffusion_4,
                   "network_1_diffusion_8" = network_1_diffusion_8,
                   "network_2_diffusion_1" = network_2_diffusion_1,
                   "network_2_diffusion_2" = network_2_diffusion_2,
                   "network_2_diffusion_4" = network_2_diffusion_4,
                   "network_2_diffusion_8" = network_2_diffusion_8,
                   "network_1_output_1" = network_1_output_1,
                   "network_1_output_2" = network_1_output_2,
                   "network_1_output_4" = network_1_output_4,
                   "network_1_output_8" = network_1_output_8,
                   "network_2_output_1" = network_2_output_1,
                   "network_2_output_2" = network_2_output_2,
                   "network_2_output_4" = network_2_output_4,
                   "network_2_output_8" = network_2_output_8,
                   "cost_matrix" = cost_matrix)    
  } else if (sum(dim(matrix_1) == dim(matrix_2)) == 2) { # Don't include the padded version if the networks are the same size
    output <- list("score" = alignment_score, 
                   "alignment" = alignment)
  } else {
    output <- list("score" = alignment_score, 
                   "alignment" = alignment, 
                   "score_with_padding" = alignment_score_padding, 
                   "alignment_with_padding" = alignment_padding)    
  }

  return(output)
}

#' @title Gini coefficient
#'
#' @description Takes a matrix and returns the Gini coefficient of each column.
#'
#' @param input A matrix where the Gini coefficient will be calculated on each column. Note that vector data must be converted to a single-column matrix.
#' 
#' @param byrow Defaults to FALSE. Set to TRUE to calculate the Gini coefficient of each row.
#'
#' @return A vector of the Gini coefficients of each column.
#' 
#' @references Gini, C. (1912). Variabilita e mutabilita. Reprinted in Memorie di metodologica statistica (Ed. Pizetti E, Salvemini, T). Rome: Libreria Eredi Virgilio Veschi.
#'
#' @examples
#' # Vectors are not supported. First convert to a single-column matrix.
#' sample_data <- runif(20, 0, 1)
#' gini(matrix(sample_data, ncol = 1))
#' 
#' # Multiple Gini coefficients can be calculated simultaneously
#' gini(matrix(sample_data, ncol = 2)) 
#' 
#' @export

gini <- function(input, byrow = FALSE)
{
  # Convert row-oriented input to columns 
  if (byrow == TRUE) {
    input <- t(input)
  }
  
  rows <- nrow(input)
  cols <- ncol(input)
  
  # Convert the input to a column-stochastic matrix (won't affect inputs already in that format)
  matrix <- sweep(input, 2, Matrix::colSums(input), FUN = "/")

  # Sort each column into ascending order
  sorted <- apply(matrix, 2, sort)

  # Turn each column into a cumulative distribution function (CDF)
  cdf <- matrix(NA, nrow = rows, ncol = cols)
  cdf[1, ] <- sorted[1, ]
  for (i in 2:rows) {
    cdf[i, ] = cdf[i - 1, ] + sorted[i, ]
  }
  cdf <- rbind(rep(0, cols), cdf)

  # The Gini coefficient is relative to a perfectly equitable distribution
  equality <- seq(from = 0, to = 1, by = 1/rows)

  # Calculate the Gini coefficient for each column (row if byrow = TRUE)
  output <- rep(NA, cols)
  for (i in 1:cols) {
    output[i] <- 2*(0.5 - pracma::trapz(equality, cdf[, i]))
  }
  
  # The output is a vector of Gini coefficients for each column (row if byrow = TRUE) of the input matrix
  return(output)
}

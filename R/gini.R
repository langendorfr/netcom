#' @title Gini coefficient
#'
#' @description Takes a matrix and returns the Gini coefficient of each row.
#'
#' @param input A row-stochastic matrix, meaning each row is a multinomial distribution.
#'
#' @return A vector of the Gini coefficient of each row.
#'
#' @examples
#' gini(matrix(runif(25, 0, 1), nrow = 5, ncol = 5))
#' netcom::gini(matrix(runif(25, 0, 1), nrow = 5, ncol = 5))
#'

gini <- function(input)
{
  # The input should be a row-stochastic matrix
  size <- nrow(input)

  sorted <- t(apply(input, 1, sort))

  cdf <- sorted
  for (i in 2:size) {
    cdf[ , i] = cdf[ , i - 1] + cdf[, i]
  }

  # The Gini coefficient is relative to a perfectly equitable distribution
  equality <- seq(from = 1/size, to = 1, by = 1/size)
  output <- rowSums(abs(sweep(cdf, 2, equality))) / size

  # The output is a vector of Gini coefficients for each row of the input matrix
  return(output)
}

#' Maximum likelihood estimation of the parameters of a matrix normal distribution
#'
#' Maximum likelihood estimation for the parameters of a matrix normal
#'    distribution \bold{X}, which is characterized by a simply separable
#'    variance-covariance structure. In the general case, which is the case
#'    considered here, two unstructured factor
#'    variance-covariance matrices determine the covariability of random
#'    matrix entries, depending on the row (one factor matrix)
#'    and the column (the other factor matrix) where two
#'    \bold{X}-entries are. In the required function, the Id1 and Id2 variables correspond to
#'    the row and column subscripts, respectively; \dQuote{value2d}
#'    indicates the observed variable.
#'
#' @param value2d from the formula value2d ~ Id1 + Id2
#' @param Id1 from the formula value2d ~ Id1 + Id2
#' @param Id2 from the formula value2d ~ Id1 + Id2
#' @param subject the replicate, also called individual
#' @param data_2d the name of the matrix data
#' @param eps the threshold in the stopping criterion for the iterative mle algorithm
#' @param maxiter the maximum number of iterations for the iterative mle algorithm
#' @param startmat the value of the second factor variance-covariance matrix used for
#'    initialization, i.e., to start the algorithm and obtain the initial estimate
#'    of the first factor variance-covariance matrix
#'
#' @section Output:
#'
#' \dQuote{Convergence}, TRUE or FALSE
#'
#' \dQuote{Iter}, will indicate the number of iterations needed for the mle algorithm to converge
#'
#' \dQuote{Xmeanhat}, the estimated mean matrix (i.e., the sample mean)
#'
#' \dQuote{First}, the row subscript, or the second column in the data file
#'
#' \dQuote{U1hat}, the estimated variance-covariance matrix for the rows
#'
#' \dQuote{Standardized.U1hat}, the standardized estimated variance-covariance matrix
#'    for the rows; the standardization is performed by dividing each entry of U1hat
#'    by entry(1, 1) of U1hat
#'
#' \dQuote{Second}, the column subscript, or the third column in the data file
#'
#' \dQuote{U2hat}, the estimated variance-covariance matrix for the columns
#'
#' \dQuote{Standardized.U2hat}, the standardized estimated variance-covariance
#'    matrix for the columns; the standardization is performed by multiplying each
#'    entry of U2hat by entry(1, 1) of U1hat
#'
#' \dQuote{Shat}, is the sample variance-covariance matrix computed from of the vectorized data matrices
#'
#' @section References:
#'
#' Dutilleul P. 1990. Apport en analyse spectrale d'un periodogramme
#' modifie et modelisation des series chronologiques avec repetitions en vue de leur comparaison
#' en frequence. D.Sc. Dissertation, Universite catholique
#' de Louvain, Departement de mathematique.
#'
#' Dutilleul P. 1999. The mle algorithm for the matrix
#' normal distribution. Journal of Statistical Computation
#' and Simulation 64: 105-123.
#'
#' @examples
#' output <- mle2d_svc(data2d$value2d, data2d$Id1, data2d$Id2, data2d$K, data_2d = data2d)
#' output
#'
#' @export

mle2d_svc <- function(value2d, Id1, Id2, subject, data_2d, eps, maxiter, startmat) {
  #Enforcing default values
  formula_2d <- value2d ~ Id1 + Id2
  data_2d$subject <- subject
  data_2d$Id1 <- Id1
  data_2d$Id2 <- Id2
  data_2d$value2d <- value2d
  data_2d <- data_2d[order(data_2d$subject, data_2d$Id1, data_2d$Id2), ]
  if (missing(eps) == TRUE) {eps <- 1e-6}
  if (missing(maxiter) == TRUE){maxiter <- 5000}
  #Quality control of data with warning for missing values, and sample size
  if (sum(is.na(data_2d)) > 0) {warning("Missing values are not accepted. Try to impute the missing values.")}
  data_2d$Id1 <- as.numeric(data_2d$Id1)
  data_2d$Id2 <- as.numeric(data_2d$Id2)
  n1 <- length(unique(data_2d$Id1))
  n2 <- length(unique(data_2d$Id2))
  K <- length(unique(data_2d$subject))
  if (missing(startmat) == TRUE){startmat <- diag(n2)}  #Value of matrix used for initialization
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (is.wholenumber(max(n1 / n2, n2 / n1))) {Kmin <- max(n1 / n2, n2 / n1) + 1} else {Kmin <- as.integer(max(n1 / n2, n2 / n1)) + 2}
  #Ensuring sufficient sample size for estimation
  if (K <= Kmin) {print("Sample size insufficient for estimation.")}
  #Setting up data for algorithm
  X <- array(data_2d$value2d, dim = c(n2, n1, K))
  Xmean <- apply(X = X, MARGIN = c(1, 2), FUN = sum) / K
  Xc <- array(0, c(n1, n2, K))
  U1int <- array(0, c(n1, n1, K))
  U2int <- array(0, c(n2, n2, K))
  #Initialization of the algorithm
  for (k in 1:K) {Xc[, , k] <- X[, , k] - Xmean}
  iter <- 0
  U2hatold <- startmat
  for (k in 1:K){U1int[, , k] <- Xc[, , k] %*% solve(U2hatold) %*% t(Xc[, , k])}
  U1hatold <- rowSums(U1int, dims = 2) / (n2 * K)
  U1int <- array(0, c(n1, n1, K))
  iter <- iter + 1
  for (k in 1:K){U2int[, , k] <- t(Xc[, , k]) %*% solve(U1hatold) %*% (Xc[, , k])}
  U2hatnew <- rowSums(U2int, dims = 2) / (n1 * K)
  U2int <- array(0, c(n2, n2, K))
  for (k in 1:K) {U1int[, , k] <- Xc[, , k] %*% solve(U2hatnew) %*% t(Xc[, , k])}
  U1hatnew <- rowSums(U1int, dims = 2) / (n2 * K)
  U1int <- array(0, c(n1, n1, K))
  #IMPORTANT: this is the MLE algorithm with iterations until convergence criterion is satisfied
  while ((norm(U1hatnew - U1hatold, "F") > eps | norm(U2hatnew-U2hatold, "F")>eps) & iter<maxiter)
  {iter=iter+1
  U1hatold <- U1hatnew
  U2hatold <- U2hatnew
  for (k in 1:K){U2int[, , k] <- t(Xc[, , k]) %*% solve(U1hatold) %*% Xc[, , k]}
  U2hatnew <- rowSums(U2int, dims = 2) / (n1 * K)
  U2int <- array(0, c(n2, n2, K))
  for (k in 1:K){U1int[, , k] <- Xc[, , k] %*% solve(U2hatnew) %*% t(Xc[, , k])}
  U1hatnew <- rowSums(U1int, dims = 2) / (n2 * K)
  U1int <- array(0, c(n1, n1, K))
  }
  if(iter == maxiter) {
    Iter <- c("Did not converge at maximum number of iterations given eps. Try to increase maxiter and/or decrease eps.")
    Convergence <- FALSE
  } else  {
    Iter <- iter
    Convergence <- TRUE}
  Shat <- U2hatnew %x% U1hatnew
  #Printing out results
  list(Call = formula_2d,
       Convergence = Convergence,
       Iter = Iter,
       Xmeanhat = Xmean,
       First = c("Id1, levels:", levels(as.factor(data_2d$Id1))),
       U1hat = U1hatnew,
       Standardized.U1hat = U1hatnew/U1hatnew[1, 1],
       Second = c("Id2, levels:", levels(as.factor(data_2d$Id2))),
       U2hat = U2hatnew,
       Standardized.U2hat = U2hatnew*(U1hatnew[1, 1]),
       Shat = Shat
  )
}

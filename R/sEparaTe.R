#' MLE and LRT functions for separable variance-covariance structures
#'
#' A package for maximum likelihood estimation (MLE) of the parameters of matrix and
#'    3rd-order tensor normal distributions with unstructured factor
#'    variance-covariance matrices (two procedures),
#'    and for unbiased modified
#'    likelihood ratio testing (LRT) of simple and double separability for
#'    variance-covariance structures (two procedures).
#'
#' @section Functions:
#'
#' mle2d_svc, for maximum likelihood estimation of the parameters of a matrix
#' normal distribution
#'
#' mle3d_svc, for maximum likelihood estimation of the parameters of a 3rd-order
#' tensor normal distribution
#'
#' lrt2d_svc, for the unbiased modified likelihood ratio test of simple
#' separability for a variance-covariance structure
#'
#' lrt3d_svc, for the unbiased modified likelihood ratio test of double
#' separability for a variance-covariance structure
#'
#' @section Data:
#'
#' data2d, a two-dimensional data set
#'
#' data3d, a three-dimensional data set
#'
#' @section References:
#'
#' Dutilleul P. 1999. The mle algorithm for the matrix
#' normal distribution. Journal of Statistical Computation
#' and Simulation 64: 105-123.
#'
#' Manceur AM, Dutilleul P. 2013. Maximum likelihood estimation for the tensor normal distribution:
#' Algorithm, minimum sample size, and empirical bias and dispersion.
#' Journal of Computational and Applied Mathematics 239: 37-49.
#'
#' Manceur AM, Dutilleul P. 2013. Unbiased modified likelihood ratio tests for simple and double
#' separability of a variance covariance structure. Statistics
#' and Probability Letters 83: 631-636.
#'
#' @docType package
#' @name sEparaTe
NULL

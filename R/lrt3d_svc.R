#' An unbiased modified likelihood ratio test for double separability of a variance-covariance structure.
#'
#' A likelihood ratio test (LRT) for double separability of a variance-covariance
#'    structure, modified to be unbiased in finite samples. The modification is a
#'    penalty-based homothetic transformation of the LRT statistic. The penalty
#'    value is optimized for a given mean model, which is left unstructured here. In
#'    the required function, the Id3, Id4 and Id5 variables correspond to the row, column
#'    and edge subscripts, respectively; \dQuote{value3d} refers to the observed
#'    variable.
#'
#' @param value3d from the formula value3d ~ Id3 + Id4 + Id5
#' @param Id3 from the formula value3d ~ Id3 + Id4 + Id5
#' @param Id4 from the formula value3d ~ Id3 + Id4 + Id5
#' @param Id5 from the formula value3d ~ Id3 + Id4 + Id5
#' @param subject the replicate, also called individual
#' @param data_3d the name of the tensor data
#' @param eps the threshold in the stopping criterion for the iterative mle
#'    algorithm (estimation)
#' @param maxiter the maximum number of iterations for the mle algorithm (estimation)
#' @param startmatU2 the value of the second factor variance-covariance
#'    matrix used for initialization
#' @param startmatU3 the value of the third factor variance-covariance matrix
#'    used for initialization, i.e., startmatU3 together with startmatU2 are used
#'    to start the mle algorithm (estimation) and obtain the initial estimate of the
#'    first factor variance-covariance matrix U1
#' @param sign.level the significance level, or rejection rate in the testing of the
#'    null hypothesis of simple separability for a variance-covariance structure, when
#'    the unbiased modified LRT is used, i.e., the critical value in the chi-square
#'    test is derived by simulations from the sampling distribution of the LRT statistic
#' @param n.simul the number of simulations used to build the sampling distribution
#'    of the LRT statistic under the null hypothesis, using the same characteristics as the
#'    i.i.d. random sample from a tensor normal distribution
#'
#' @section Output:
#'
#' \dQuote{Convergence}, TRUE or FALSE
#'
#' \dQuote{chi.df}, the theoretical number of degrees of freedom of the asymptotic
#'    chi-square distribution that would apply to the unmodified LRT statistic for double separability
#'    of a variance-covariance structure
#' \dQuote{Lambda}, the observed value of the unmodified LRT statistic
#'
#' \dQuote{critical.value}, the critical value at the specified significance
#'    level for the chi-square distribution with \dQuote{chi.df} degrees of freedom
#'
#' \dQuote{Decision.lambda}, the decision (acceptance/rejection) regarding the null
#'    hypothesis of double separability, made using the theoretical (biased unmodified) LRT
#'
#' \dQuote{Simulation.critical.value}, the critical value at the specified significance
#'    level that is derived from the sampling distribution of the unbiased modified LRT statistic
#'
#' \dQuote{Decision.lambda.simulation}, the decision (acceptance/rejection) regarding
#'    the null hypothesis of double separability, made using the unbiased modified LRT
#'
#' \dQuote{Penalty}, the optimized penalty value used in the homothetic transformation
#'    between the biased unmodified and unbiased modified LRT statistics
#'
#' \dQuote{U1hat}, the estimated variance-covariance matrix for the rows
#'
#' \dQuote{Standardized_U1hat}, the standardized estimated variance-covariance matrix
#'    for the rows; the standardization is performed by dividing each entry of U1hat by
#'    entry(1, 1) of U1hat
#'
#' \dQuote{U2hat}, the estimated variance-covariance matrix for the columns
#'
#' \dQuote{Standardized_U2hat}, the standardized estimated variance-covariance matrix
#'    for the columns; the standardization is performed by multiplying each entry of
#'    U2hat by entry(1, 1) of U1hat
#'
#' \dQuote{U3hat}, the estimated variance-covariance matrix for the edges
#'
#' \dQuote{Shat}, the sample variance-covariance matrix computed from the
#'    vectorized data tensors
#'
#' @section References:
#'
#' Manceur AM, Dutilleul P. 2013. Unbiased modified likelihood ratio tests for simple and double
#' separability of a variance-covariance structure. Statistics
#' and Probability Letters 83: 631-636.
#'
#' @examples
#' output <- lrt3d_svc(data3d$value3d, data3d$Id3, data3d$Id4, data3d$Id5,
#'             data3d$K, data_3d = data3d, n.simul = 100)
#' output
#'
#' @export

lrt3d_svc <- function(value3d, Id3, Id4, Id5 , subject, data_3d, eps, maxiter, startmatU2, startmatU3, sign.level, n.simul)
{#Enforcing default values
  formula_3d <- value3d ~ Id3 + Id4 + Id5
  data_3d$subject <- subject
  data_3d$Id3 <- Id3
  data_3d$Id4 <- Id4
  data_3d$Id5 <- Id5
  data_3d$value3d <- value3d
  data_3d <- data_3d[order(data_3d$subject, data_3d$Id3, data_3d$Id4, data_3d$Id5), ]
  if (missing(eps) == TRUE) {eps <- 1e-6}
  if (missing(maxiter) == TRUE) {maxiter <- 100}
  if (missing(sign.level) == TRUE) {sign.level <- 0.95}
  if (missing(n.simul) == TRUE) {n.simul <- 8000}
  #Quality control of data with warning for missing values, and sample size
  if (sum(is.na(data_3d)) > 0) {warning("Missing values are not accepted. Try to impute the missing values.")}
  data_3d$Id3 <- as.numeric(data_3d$Id3)
  data_3d$Id4 <- as.numeric(data_3d$Id4)
  data_3d$Id5 <- as.numeric(data_3d$Id5)
  n1 <- length(unique(data_3d$Id3))
  n2 <- length(unique(data_3d$Id4))
  n3 <- length(unique(data_3d$Id5))
  K <- length(unique(data_3d$subject))
  if (missing(startmatU2) == TRUE) {startmatU2 <- diag(n2)}  #Value of matrix used for initialization
  if (missing(startmatU3) == TRUE){startmatU3 <- diag(n3)}  #Value of matrix used for initialization
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if (is.wholenumber(max(n1 / (n2 * n3), n2 / (n1 * n3), n3 / (n1 * n2)))) {Kmin = max(n1 / (n2 * n3), n2 / (n1 * n3), n3 / (n1 * n2)) + 1} else {Kmin = as.integer(max(n1 / (n2 * n3), n2 / (n1 * n3), n3 / (n1 * n2))) + 2}
  if (K <= Kmin) {print("Sample size insufficient for estimation.")} #Ensuring sufficient sample size for estimation
  if (K < (n1 * n2 * n3) + 1) {print ("Sample size insufficient for LRT of separability")} #Ensuring sufficient sample size for estimation
  chi.df <- (n1 * n2 * n3 * (((n1 * n2 * n3) + 1) / 2)) - (n1 * (((n1) + 1) / 2)) - (n2 * (((n2) + 1) / 2)) - (n3 * (((n3) + 1) / 2)) + 1
  #Setting up data for algorithm
  X <- array(data_3d$value3d, dim = c(n3, n2, n1, K))
  Xmean <- apply(X, MARGIN = c(1, 2, 3), sum) / K
  Xc <- array(0, dim = c(n3, n2, n1, K))
  U1int <- array(0, dim = c(n1, n1, K))
  U2int <- array(0, dim = c(n2, n2, K))
  U3int <- array(0, dim = c(n3, n3, K))
  for (k in 1:K) {Xc[, , , k] <- X[, , , k] - Xmean}
  XU1 <- array(aperm(Xc, perm = c(1, 2, 3, 4)), dim = c(n1, n3 * n2, K))
  XU2 <- array(aperm(Xc, perm = c(2, 1, 3, 4)), dim = c(n2, n3 * n1, K))
  XU3 <- array(aperm(Xc, perm = c(3, 1, 2, 4)), dim = c(n3, n2 * n1, K))
  #Initialization of the algorithm
  iter <- 0
  U3hatold <- startmatU3
  U2hatold <- startmatU2
  tt1 <- U3hatold %x% U2hatold
  for (k in 1:K) {U1int[ , , k] <- XU1[, , k] %*% solve(tt1) %*% aperm(XU1[, , k], perm = c(2, 1))}
  U1hatold <- apply(U1int, MARGIN = c(1, 2), sum) / (n2 * n3 * K)
  iter <- iter + 1
  tt2 <- U3hatold %x% U1hatold
  for (k in 1:K){U2int[, , k] <- XU2[, , k] %*% solve(tt2) %*% aperm(XU2[, , k], perm = c(2, 1))}
  U2hatnew <- apply(U2int, MARGIN = c(1, 2), sum) / (n3 * n1 * K)
  tt3 <- U2hatnew %x% U1hatold
  for (k in 1:K){U3int[, , k] <- XU3[, , k] %*% solve(tt3) %*% aperm(XU3[, , k], perm = c(2, 1))}
  U3hatnew <- apply(U3int, MARGIN = c(1, 2), sum) / (n2 * n1 * K)
  tt1 <- U3hatnew %x% U2hatnew
  for (k in 1:K) {U1int[, , k] <- XU1[, , k] %*% solve(tt1) %*% aperm(XU1[, , k], perm = c(2, 1))}
  U1hatnew <- apply(U1int, MARGIN = c(1, 2), sum) / (n2 * n3 * K)
  #IMPORTANT: this is the MLE algorithm with iterations until convergence criterion is satisfied
  while ((norm(U1hatnew - U1hatold, "F") > eps | norm(U2hatnew - U2hatold, "F") > eps | norm(U3hatnew - U3hatold, "F") > eps)  & iter < maxiter)
  {iter = iter + 1
  U1hatold <- U1hatnew
  U2hatold <- U2hatnew
  U3hatold <- U3hatnew
  tt2 <- U3hatold %x% U1hatold
  for (k in 1:K){U2int[, , k] <- XU2[, , k] %*% solve(tt2) %*% aperm(XU2[, , k], perm=c(2, 1))}
  U2hatnew <- apply(U2int, MARGIN = c(1, 2), sum) / (n3 * n1 * K)
  tt3 <- U2hatnew %x% U1hatold
  for (k in 1:K){U3int[, , k] <- XU3[, , k] %*% solve(tt3) %*% aperm(XU3[, , k], perm = c(2, 1))}
  U3hatnew <- apply(U3int, MARGIN = c(1, 2), sum) / (n2 * n1 * K)
  tt1 <- U3hatnew %x% U2hatnew
  for (k in 1:K) {U1int[, , k] <- XU1[, , k] %*% solve(tt1) %*% aperm(XU1[, , k], perm = c(2, 1))}
  U1hatnew <- apply(U1int, MARGIN = c(1, 2), sum) / (n2 * n3 * K)
  }
  if(iter == maxiter) {
    Iter <- c("Did not converge at maximum number of iterations given eps. Try to increase maxiter and/or decrease eps.")
    Convergence <- FALSE
  } else  {
    Iter <- iter
    Convergence <- TRUE}
    Shat <- U3hatnew %x% U2hatnew %x% U1hatnew
  #ML estimation under vector normal model
  X.un <- t(apply(X, MARGIN = c(4), function(x) as.vector(x)))
  Xmean.un <- apply(X.un, 2, sum) / K
  S.un.all <- array(0, c(n1 * n2 * n3, n1 * n2 * n3, K))
  for (k in 1:K) {S.un.all[, , k] <- as.matrix(X.un[k, ] - Xmean.un) %*% t(as.matrix(X.un[k, ] - Xmean.un))}
  S.un <- rowSums(S.un.all, dims = 2) / K
  #Lambda
  Lambda <- (K - 1) * ((n1 * n2 * log10(det(U3hatnew))) + (n1 * n3 * log10(det(U2hatnew))) + (n2 * n3 * log10(det(U1hatnew))) - log10(det(S.un)))
  U1out <- U1hatnew
  U2out <- U2hatnew
  U3out <- U3hatnew
  theoretical.quantiles <- as.numeric(stats::quantile(Lambda, probs = c(.01, .05, .1, .25, .5, .75, .9, .95, .99)))
  #Lambda star
  trueU1 <- diag(n1)
  trueU2 <- diag(n2)
  trueU3 <- diag(n3)
  trueS <- trueU3 %x% trueU2 %x% trueU1
  simulS <- t(chol(trueS))
  Lambda.simul <- rep(0, n.simul)
  maxiter.simul <- 10000
  for (i in 1:n.simul){
    Xsimul.all <- list()
    for (k in 1:K){
      Xsimul.loop <- simulS %*% stats::rnorm(n1 * n2 * n3) #This is just the vector of z values
      Xsimul.all[[k]] <- Xsimul.loop
    }
    Xsimul <- array(unlist(Xsimul.all), dim = c(n3, n2, n1, K))
    #ML estimation under vector normal model
    X.un <- t(apply(Xsimul, MARGIN = c(4), function(x) as.vector(x)))
    Xmean.un <- apply(X.un, 2, sum) / K
    S.un.all <- array(0, c(n1*n2*n3, n1*n2*n3, K))
    for (k in 1:K) {S.un.all[, , k] <- as.matrix(X.un[k, ] - Xmean.un) %*% t(as.matrix(X.un[k, ] - Xmean.un))}
    S.un.simul <- rowSums(S.un.all, dims = 2) / K
    #Setting up data for algorithm
    Xmean <- apply(Xsimul, MARGIN = c(1, 2, 3), sum) / K
    Xc <- array(0, dim = c(n3, n2, n1, K))
    U1int <- array(0, dim = c(n1, n1, K))
    U2int <- array(0, dim = c(n2, n2, K))
    U3int <- array(0, dim = c(n3, n3, K))
    for (k in 1:K) {Xc[, , , k] <- Xsimul[, , , k] - Xmean}
    XU1 <- array(aperm(Xc, perm = c(1, 2, 3, 4)), dim = c(n1, n3 * n2, K))
    XU2 <- array(aperm(Xc, perm = c(2, 1, 3, 4)), dim = c(n2, n3 * n1, K))
    XU3 <- array(aperm(Xc, perm = c(3, 1, 2, 4)), dim = c(n3, n2 * n1, K))
    #Initialization of the algorithm
    iter <- 0
    U3hatold <- startmatU3
    U2hatold <- startmatU2
    tt1 <- U3hatold %x% U2hatold
    for (k in 1:K) {U1int[, , k] <- XU1[, , k] %*% solve(tt1) %*% aperm(XU1[, , k], perm = c(2, 1))}
    U1hatold <- apply(U1int, MARGIN = c(1,2), sum) / (n2 * n3 * K)
    iter <- iter + 1
    tt2 <- U3hatold %x% U1hatold
    for (k in 1:K){U2int[, , k] <- XU2[, , k] %*% solve(tt2) %*% aperm(XU2[, , k], perm = c(2, 1))}
    U2hatnew <- apply(U2int, MARGIN = c(1,2), sum) / (n3 * n1 * K)
    tt3 <- U2hatnew %x% U1hatold
    for (k in 1:K){U3int[, , k] <- XU3[, , k] %*% solve(tt3) %*% aperm(XU3[, , k], perm = c(2, 1))}
    U3hatnew <- apply(U3int, MARGIN = c(1, 2), sum) / (n2 * n1 * K)
    tt1 <- U3hatnew %x% U2hatnew
    for (k in 1:K) {U1int[, , k] <- XU1[, , k] %*% solve(tt1) %*% aperm(XU1[, , k], perm = c(2, 1))}
    U1hatnew <- apply(U1int, MARGIN = c(1, 2), sum) / (n2 * n3 * K)
    #IMPORTANT: this is the MLE algorithm with iterations until convergence criterion is satisfied
    while ((norm(U1hatnew - U1hatold, "F") > eps | norm(U2hatnew - U2hatold, "F") > eps | norm(U3hatnew - U3hatold, "F") > eps) & iter < maxiter.simul)
    {iter <- iter + 1
    U1hatold <- U1hatnew
    U2hatold <- U2hatnew
    U3hatold <- U3hatnew
    tt2 <- U3hatold %x% U1hatold
    for (k in 1:K) {U2int[, , k] <- XU2[, , k] %*% solve(tt2) %*% aperm(XU2[, , k], perm = c(2, 1))}
    U2hatnew <- apply(U2int, MARGIN = c(1, 2), sum) / (n3 * n1 * K)
    tt3 <- U2hatnew %x% U1hatold
    for (k in 1:K) {U3int[, , k] <- XU3[, , k] %*% solve(tt3) %*% aperm(XU3[, , k], perm = c(2, 1))}
    U3hatnew <- apply(U3int, MARGIN = c(1, 2), sum) / (n2 * n1 * K)
    tt1 <- U3hatnew %x% U2hatnew
    for (k in 1:K) {U1int[, , k] <- XU1[, , k] %*% solve(tt1) %*% aperm(XU1[, , k], perm = c(2, 1))}
    U1hatnew <- apply(U1int, MARGIN = c(1, 2), sum) / (n2 * n3 * K)
    }
    U1hatnew.simul <- U1hatnew
    U2hatnew.simul <- U2hatnew
    U3hatnew.simul <- U3hatnew
    Lambda.simul[i] <- (K - 1) * ((n1 * n2 * log10(det(U3hatnew.simul))) + (n1 * n3 * log10(det(U2hatnew.simul))) + (n2 * n3 * log10(det(U1hatnew.simul))) - log10(det(S.un.simul)))
  }
  modified.critical.value <- as.numeric(stats::quantile(Lambda.simul, sign.level))
  empirical.quantiles <- as.numeric(stats::quantile(Lambda.simul, prob = c(.01, .05, .1, .25, .5, .75, .9, .95, .99)))
  #Decision
  critical.value <- stats::qchisq(sign.level, chi.df)
  Decision.lambda <- c()
  if (Lambda > critical.value){
    Decision.lambda <- c("Reject null hypothesis of separability.")
  } else {
    Decision.lambda <- c("Fail to reject null hypothesis of separability.")
  }
  Decision.lambda.modified <- c()
  if (Lambda > modified.critical.value){
    Decision.lambda.modified <- c("Reject null hypothesis of separabiltiy.")
  } else {
    Decision.lambda.modified <- c("Fail to reject null hypothesis of separability.")
  }
  #Printing out results
  list(Convergence = Convergence,
       chi.df = chi.df,
       Lambda = Lambda,
       critical.value = critical.value,
       Decision.lambda = Decision.lambda,
       Simulation.critical.value = modified.critical.value,
       Penalty = as.numeric(-(length(unique(data_3d$subject))) * (stats::coef(stats::lm(empirical.quantiles ~ theoretical.quantiles - 1))) + (length(unique(data_3d$subject)))),
       Decision.lambda.simulation = Decision.lambda.modified,
       U1hat = U1out,
       Standardized.U1hat = U1out / U1out[1, 1],
       U2hat = U2out,
       Standardized.U2hat = U2out * (U1out[1, 1]),
       U3hat = U3out,
       Shat = Shat
  )
}


#' Unbiased modified likelihood ratio test for simple separability of a variance-covariance matrix.
#'
#' A likelihood ratio test (LRT) for simple separability of a
#' variance-covariance matrix, modified to be unbiased in finite samples. The
#' modification is a penalty-based homothetic transformation of the LRT
#' statistic. The penalty value is optimized for a given mean model, which is
#' left unstructured here. In the required function, the Id1 and Id2 variables
#' correspond to the row and column subscripts, respectively; \dQuote{value2d} refers to
#' the observed variable.
#'
#' @param value2d from the formula value2d ~ Id1 + Id2
#' @param Id1 from the formula value2d ~ Id1 + Id2
#' @param Id2 from the formula value2d ~ Id1 + Id2
#' @param subject the replicate, also called the subject or individual, the
#'    first column in the matrix (2d) data file
#' @param data_2d the name of the matrix data
#' @param eps the threshold in the stopping criterion for the iterative mle
#'    algorithm (estimation)
#' @param maxiter the maximum number of iterations for the mle algorithm (estimation)
#' @param startmat the value of the second factor variance-covariance matrix
#'    used for initialization, i.e., to start the mle algorithm (estimation) and
#'    obtain the initial estimate of the first factor variance-covariance matrix
#' @param sign.level the significance level, or rejection rate in the testing of
#'    the null hypothesis of simple separability for a variance-covariance
#'    structure, when the unbiased modified LRT is used, i.e., the critical value
#'    in the chi-square test is derived by simulations from the sampling distribution
#'    of the LRT statistic
#' @param n.simul the number of simulations used to build the sampling distribution
#'    of the LRT statistic under the null hypothesis, using the same characteristics as the
#'    i.i.d. random sample from a matrix normal distribution
#'
#' @section Output:
#'
#' \dQuote{Convergence}, TRUE or FALSE
#'
#' \dQuote{chi.df}, the theoretical number of degrees of freedom of the
#'    asymptotic chi-square distribution that would apply to the unmodified
#'    LRT statistic for simple
#'    separability of a variance-covariance structure
#'
#' \dQuote{Lambda}, the observed value of the unmodified LRT statistic
#'
#' \dQuote{critical.value}, the critical value at the specified significance
#'    level for the chi-square distribution with \dQuote{chi.df} degrees of freedom
#'
#' \dQuote{Decision.lambda} will indicate whether or not the null hypothesis of
#'    separability was rejected, based on the theoretical LRT statistic
#'
#' \dQuote{Simulation.critical.value}, the critical value at the specified
#'    significance level that is derived from the sampling distribution of the
#'    unbiased modified LRT statistic
#'
#' \dQuote{Decision.lambda.simulation}, the decision (acceptance/rejection)
#'    regarding the null hypothesis of simple separability, made using the
#'    theoretical (biased unmodified) LRT
#'
#' \dQuote{Penalty}, the optimized penalty value used in the homothetic
#'    transformation between the biased unmodified and unbiased modified LRT statistics
#'
#' \dQuote{U1hat}, the estimated variance-covariance matrix for the rows
#'
#' \dQuote{Standardized_U1hat}, the standardized estimated variance-covariance
#'    matrix for the rows; the standardization is performed by dividing each
#'    entry of U1hat by entry(1, 1) of U1hat
#'
#' \dQuote{U2hat}, the estimated variance-covariance matrix for the columns
#'
#' \dQuote{Standardized_U2hat}, the standardized estimated variance-covariance
#'    matrix for the columns; the standardization is performed by multiplying
#'    each entry of U2hat by entry(1, 1) of U1hat
#'
#' \dQuote{Shat}, the sample variance-covariance matrix computed from the
#'    vectorized data matrices

#' @section References:
#'
#' Manceur AM, Dutilleul P. 2013. Unbiased modified likelihood ratio tests for simple and double
#' separability of a variance-covariance structure. Statistics
#' and Probability Letters 83: 631-636.
#'
#' @examples
#' output <- lrt2d_svc(data2d$value2d, data2d$Id1, data2d$Id2,
#'             data2d$K, data_2d = data2d, n.simul = 100)
#' output
#'
#' @export

lrt2d_svc <- function(value2d, Id1, Id2, subject, data_2d, eps, maxiter, startmat, sign.level, n.simul) {
#Enforcing default values
  formula_2d <- value2d ~ Id1 + Id2
  data_2d$subject <- subject
  data_2d$Id1 <- Id1
  data_2d$Id2 <- Id2
  data_2d$value2d <- value2d
  data_2d <- data_2d[order(data_2d$subject, data_2d$Id1, data_2d$Id2), ]
  if (missing(eps) == TRUE){eps <- 1e-6}
  if (missing(maxiter) == TRUE){maxiter <- 5000}
  if (missing(sign.level) == TRUE){sign.level <- 0.95}
  if (missing(n.simul) == TRUE) {n.simul <- 8000}
  #Quality control of data with warning for missing values, and sample size
  mf <- stats::model.frame(formula = formula_2d, data = data_2d)
  if (sum(is.na(data_2d)) > 0){warning("Missing values are not accepted. Try to impute the missing values.")}
  data_2d$Id1 <- as.numeric(data_2d$Id1)
  data_2d$Id2 <- as.numeric(data_2d$Id1)
  n1 <- length(unique(data_2d$Id1))
  n2 <- length(unique(data_2d$Id2))
  K <- length(unique(data_2d$subject))
  if (missing(startmat) == TRUE){startmat <- diag(n2)}  #Value of matrix used for initialization
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if (is.wholenumber(max(n1 / n2, n2 / n1))) {Kmin <- max(n1 / n2, n2 / n1) + 1} else {Kmin <- as.integer(max(n1 / n2, n2 / n1)) + 2}
  if (K <= Kmin) {print("Sample size insufficient for estimation.")} #Ensuring sufficient sample size for estimation
  if (K < (n1 * n2) + 1) {print ("Sample size insufficient for LRT of separability")} #Ensuring sufficient sample size for estimation
  chi.df <- (n1 * n2 * (((n1 * n2) + 1) / 2)) - (n1 * (((n1) + 1) / 2)) - (n2 * (((n2) + 1) / 2)) + 1
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
  for (k in 1:K){U2int[, , k] <- t(Xc[, , k]) %*% solve(U1hatold) %*% Xc[, , k]}
  U2hatnew <- rowSums(U2int, dims = 2) / (n1 * K)
  U2int <- array(0, c(n2, n2, K))
  for (k in 1:K){U1int[, , k] <- Xc[, , k] %*% solve(U2hatnew) %*% t(Xc[, , k])}
  U1hatnew <- rowSums(U1int, dims = 2) / (n2 * K)
  U1int <- array(0, c(n1, n1, K))
  #IMPORTANT: this is the MLE algorithm with iterations until convergence criterion is satisfied
  while ((norm(U1hatnew - U1hatold, "F") > eps | norm(U2hatnew - U2hatold, "F") > eps) & iter < maxiter)
  {iter <- iter + 1
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
  } else {
    Iter <- iter
    Convergence <- TRUE}
  #ML estimation under vector normal model
  X.un <- matrix(0, K, n1*n2)
  for (k in 1:K) X.un[k, ] <- as.vector(X[, , k])
  Xmean.un <- apply(X.un, 2, sum) / K
  S.un.all <- array(0, c(n1 * n2, n1 * n2, K))
  for (k in 1:K){S.un.all[, , k] <- as.matrix(X.un[k, ] - Xmean.un) %*% t(as.matrix(X.un[k, ] - Xmean.un))}
  S.un <- rowSums(S.un.all, dims = 2) / K
  #Lambda
  Lambda <- (K - 1) * ((n2 * log10(det(U1hatnew))) + (n1 * log10(det(U2hatnew))) - log10(det(S.un)))
  U1out <- U1hatnew
  U2out <- U2hatnew
  #Lambda star
  mat.sqrt <- function(A)
  {
    ei <- eigen(A)
    d <- ei$values
    d <- (d + abs(d)) / 2
    d2 <- sqrt(d)
    ans <- ei$vectors %*% diag(d2) %*% t(ei$vectors)
  }
  trueU1 <- diag(n1)
  trueU2 <- diag(n2)
  Lambda.simul <- rep(0, n.simul)
  for (i in 1:n.simul) {
    Xsimul.all <- list()
    for (k in 1:K){
      Xsimul <- ((mat.sqrt(trueU1)) %*% (matrix(stats::rnorm(n1 * n2), n1, n2)) %*% (mat.sqrt(trueU2)))
      infoX <- matrix(c(rep(1:n1, each = n2), rep(1:n2, n1), rep(k, n1 * n2)), nrow = n1 * n2 , ncol = 3)
      Xsimul <- matrix(t(Xsimul), nrow = n1 * n2, ncol = 1)
      Xsimul.all[[k]] <- cbind(Xsimul, infoX)
    }
    Xsimul.all <- do.call(rbind.data.frame, Xsimul.all)
    names(Xsimul.all) <- c("value", "Id1", "Id2", "K")
    dataXsimul.all <- split(Xsimul.all, Xsimul.all[,dim(Xsimul.all)[2]])
    Xsimul.all <- lapply(dataXsimul.all, function(x) t(matrix(x$value, nrow = n2, ncol = n1)))
        #UNSTRUCTURED
    X.un.simul <- lapply(Xsimul.all, function(x) as.vector(x))
    Xmean.un.simul <- Reduce("+", X.un.simul) / K
    S.un.all.simul <- array(0, c(n1*n2, n1*n2, K))
    for (k in 1:K){S.un.all.simul[, , k] <- as.matrix(X.un.simul[[k]] - Xmean.un.simul) %*% t(as.matrix(X.un.simul[[k]] - Xmean.un.simul))}
    S.un.simul <- rowSums(S.un.all.simul, dims = 2) / K
    #SEPARABLE
    X <- lapply(dataXsimul.all, function(x) t(matrix(x$value, nrow = n2, ncol = n1)))
    Xmean <- Reduce("+", X) / K
    Xc <- array(0, c(n1, n2, K))
    U1int <- array(0, c(n1, n1, K))
    U2int <- array(0, c(n2, n2, K))
        #Initialization of the algorithm
    for (k in 1:K) {Xc[, , k] <- X[[k]] - Xmean}
    iter <- 0
    U2hatold <- startmat
    for (k in 1:K) {U1int[, , k] <- Xc[, , k] %*% solve(U2hatold) %*% t(Xc[, , k])}
    U1hatold <- rowSums(U1int, dims = 2)/ (n2 * K)
    U1int <- array(0, c(n1, n1, K))
    iter <- iter + 1
    for (k in 1:K){U2int[, , k] <- t(Xc[, , k]) %*% solve(U1hatold) %*% Xc[, , k]}
    U2hatnew <- rowSums(U2int, dims = 2) / (n1 * K)
    U2int <- array(0, c(n2, n2, K))
    for (k in 1:K){U1int[, , k] <- Xc[, , k] %*% solve(U2hatnew) %*% t(Xc[, , k])}
    U1hatnew <- rowSums(U1int, dims = 2) / (n2 * K)
    U1int <- array(0, c(n1, n1, K))
    #IMPORTANT: this is the MLE algorithm with iterations until convergence criterion is satisfied
    while ((norm(U1hatnew - U1hatold, "F") > eps | norm(U2hatnew - U2hatold, "F") > eps) & iter < maxiter)
    {iter <- iter + 1
    U1hatold <- U1hatnew
    U2hatold <- U2hatnew
    for (k in 1:K) {U2int[, , k] <- t(Xc[, , k]) %*% solve(U1hatold) %*% Xc[, , k]}
    U2hatnew <- rowSums(U2int, dims = 2) / (n1 * K)
    #U2int<-c()
    U1int <- array(0, c(n1, n1, K))
    for (k in 1:K) {U1int[, , k] <- Xc[, , k] %*% solve(U2hatnew) %*% t(Xc[, , k])}
    U1hatnew <- rowSums(U1int, dims = 2) / (n2 * K)
    U1int <- array(0, c(n1, n1, K))
    }
    Lambda.simul[i] <- (K-1) * ((n2 * log10(det(U1hatnew))) + (n1*log10(det(U2hatnew))) - log10(det(S.un.simul)))
  }
  modified.critical.value <- as.numeric(stats::quantile(Lambda.simul, sign.level))
  empirical.quantiles <- as.numeric(stats::quantile(Lambda.simul, prob = c(.01, .05, .1, .25, .5, .75, .9, .95, .99)))
  theoretical.quantiles <- as.numeric(stats::quantile(Lambda, prob = c(.01, .05, .1, .25, .5, .75, .9, .95, .99)))
  #Decision
  critical.value <- stats::qchisq(sign.level, chi.df)
  Decision.lambda <- c()
  if (Lambda > critical.value){
    Decision.lambda <- c("Reject null hypothesis of separabiltiy.")
  } else {
    Decision.lambda <- c("Fail to reject null hypothesis of separability.")
  }
  Decision.lambda.modified <- c()
  if (Lambda > modified.critical.value){
    Decision.lambda.modified<-c("Reject null hypothesis of separabiltiy.")
  } else {
    Decision.lambda.modified<-c("Fail to reject null hypothesis of separability.")
  }
  list(Convergence = Convergence,
       chi.df = chi.df,
       Lambda = Lambda,
       critical.value = critical.value,
       Decision.lambda = Decision.lambda,
       Simulation.critical.value = modified.critical.value,
       Decision.lambda.simulation = Decision.lambda.modified,
       Penalty = as.numeric(-(length(unique(subject))) * (stats::coef(stats::lm(empirical.quantiles ~ theoretical.quantiles - 1))) + (length(unique(subject)))),
       U1hat = U1out,
       Standardized.U1hat = U1out / U1out[1, 1],
       U2hat = U2out,
       Standardized.U2hat = U2out * U1out[1, 1],
       Shat = (U2out * U1out[1, 1]) %x% (U1out / U1out[1, 1])
  )
}


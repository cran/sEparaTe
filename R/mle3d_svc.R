#' Maximum likelihood estimation of the parameters of a 3rd-order tensor normal distribution
#'
#' Maximum likelihood estimation for the parameters of a
#' 3rd-order tensor normal distribution \strong{X}, which is characterized by
#' a doubly separable variance-covariance structure. In the general
#' case, which is the case considered here, three unstructured
#' factor variance-covariance matrices determine the covariability
#' of random tensor entries, depending on the row (one factor
#' matrix), the column (another factor matrix) and the edge (remaining
#' factor matrix) where two \strong{X}-entries are. In the required function, the
#' Id3, Id4 and Id5 variables correspond to the row, column and edge
#' subscripts, and are the second, third and fourth columns in the tensor
#' (3d) data file, respectively; \dQuote{value3d} indicates the observed
#' variable, and is the fifth column in the tensor data file.
#'
#' @param formula value3d~Id3+Id4+Id5
#' @param rep the replicate, also called subject or individual, the first column in the tensor (3d) data file
#' @param data the name of the tensor data
#' @param eps the threshold in the stopping criterion for the iterative mle algorithm
#' @param maxiter the maximum number of iterations for the iterative mle algorithm
#' @param startmatU2 the value of the second factor variance covariance matrix used for initialization
#' @param startmatU3 the value of the third factor variance covariance matrix used for initialization,
#'    i.e., startmatU3 together with startmatU2 are used to start the algorithm and obtain the initial estimate of
#'    the first factor variance covariance matrix U1
#'
#' @section Output:
#'
#' \dQuote{Convergence}, TRUE or FALSE
#'
#' \dQuote{Iter}, the number of iterations needed for the mle algorithm to converge
#'
#' \dQuote{Xmeanhat}, the estimated mean tensor (i.e., the sample mean)
#'
#' \dQuote{First}, the row subscript, or the second column in the data file
#'
#' \dQuote{U1hat}, the estimated variance-covariance matrix for the rows
#'
#' \dQuote{Standardized.U1hat}, the standardized estimated variance-covariance matrix for the rows; the standardization is performed by dividing each entry of U1hat by entry(1, 1) of U1hat
#'
#' \dQuote{Second}, the column subscript, or the third column in the data file
#'
#' \dQuote{U2hat}, the estimated variance-covariance matrix for the columns
#'
#' \dQuote{Standardized.U2hat}, the standardized estimated variance-covariance matrix for the columns; the standardization is performed by multiplying each entry of U2hat by entry(1, 1) of U1hat
#'
#' \dQuote{Third}, the edge subscript, or the fourth column in the data file
#'
#' \dQuote{U3hat}, the estimated variance-covariance matrix for the edges
#'
#' \dQuote{Shat}, the sample variance-covariance matrix computed from the vectorized data tensors
#'
#' @section Reference:
#'
#' Manceur AM, Dutilleul P. 2013. Maximum likelihood estimation for the tensor normal distribution:
#' Algorithm, minimum sample size, and empirical bias and dispersion.
#' Journal of Computational and Applied Mathematics 239: 37-49.
#'
#' @examples
#' output <- mle3d_svc(value3d~Id3+Id4+Id5, rep = "K", data = data3d)
#' output
#'
#' @export

mle3d_svc <- function(formula, rep, data=list(), eps, maxiter, startmatU2, startmatU3)
{#Enforcing default values
  formula<-paste(deparse(formula), eval(rep), sep="+")
  formula<-stats::as.formula(formula)
  if (missing(eps)==TRUE){eps=1e-6}
  if (missing(maxiter)==TRUE){maxiter=100}
  #Quality control of data with warning for missing values, and sample size
  mf <- stats::model.frame(formula=formula, data=data)
  if (sum(is.na(mf))>0){warning("Missing values are not accepted. Try to impute the missing values.")}
  colnames(mf)[1]<-"value"
  rpl.mf<-c(dim(mf)[2])
  mf <-mf[order(mf[, rpl.mf]),]
  mf[,2]<-as.numeric(mf[,2])
  mf[,3]<-as.numeric(mf[,3])
  mf[,4]<-as.numeric(mf[,4])
  n1<-length(unique(mf[,2]))
  n2<-length(unique(mf[,3]))
  n3<-length(unique(mf[,4]))
  K<-length(unique(mf[,rpl.mf]))
  if (missing(startmatU2)==TRUE){startmatU2<-diag(n2)}  #Value of matrix used for initialization
  if (missing(startmatU3)==TRUE){startmatU3<-diag(n3)}  #Value of matrix used for initialization
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (is.wholenumber(max(n1/(n2*n3), n2/(n1*n3), n3/(n1*n2)))) {Kmin=max(n1/(n2*n3), n2/(n1*n3), n3/(n1*n2))+1} else {Kmin=as.integer(max(n1/(n2*n3), n2/(n1*n3), n3/(n1*n2)))+2}
  if (K<=Kmin) {print("Sample size insufficient for estimation.")} #Ensuring sufficient sample size for estimation
  #Setting up data for algorithm
  X<-array(mf$value, dim=c(n3,n2,n1,K))
  Xmean<-apply(X, MARGIN=c(1,2,3), sum)/K
  Xc<-array(0,dim=c(n3,n2,n1,K))
  U1int<-array(0, dim=c(n1,n1,K))
  U2int<-array(0, dim=c(n2,n2,K))
  U3int<-array(0, dim=c(n3,n3,K))
  for (k in 1:K) {Xc[,,,k]<-X[,,,k]-Xmean}
  XU1=array(aperm(Xc, perm=c(1,2,3,4)), dim=c(n1, n3*n2, K))
  XU2=array(aperm(Xc, perm=c(2,1,3,4)), dim=c(n2, n3*n1, K))
  XU3=array(aperm(Xc, perm=c(3,1,2,4)), dim=c(n3, n2*n1, K))
  U3hatold=startmatU2
  U2hatold=startmatU3
  iter=0
  tt1=U3hatold%x%U2hatold
  for (k in 1:K){U1int[,,k]=XU1[,,k]%*%solve(tt1)%*%aperm(XU1[,,k], perm=c(2,1))}
  U1hatold<-apply(U1int, MARGIN=c(1,2), sum)/(n2*n3*K)
  iter=iter+1
  tt2=U3hatold%x%U1hatold
  for (k in 1:K){U2int[,,k]=XU2[,,k]%*%solve(tt2)%*%aperm(XU2[,,k], perm=c(2,1))}
  U2hatnew<-apply(U2int, MARGIN=c(1,2), sum)/(n3*n1*K)
  tt3=U2hatnew%x%U1hatold
  for (k in 1:K){U3int[,,k]=XU3[,,k]%*%solve(tt3)%*%aperm(XU3[,,k], perm=c(2,1))}
  U3hatnew<-apply(U3int, MARGIN=c(1,2), sum)/(n2*n1*K)
  tt1=U3hatnew%x%U2hatnew
  for (k in 1:K){U1int[,,k]=XU1[,,k]%*%solve(tt1)%*%aperm(XU1[,,k], perm=c(2,1))}
  U1hatnew<-apply(U1int, MARGIN=c(1,2), sum)/(n2*n3*K)
  #IMPORTANT: this is the MLE algorithm with iterations until convergence criterion is satisfied
  while ((norm(U1hatnew-U1hatold, "F")>eps | norm(U2hatnew-U2hatold, "F")>eps | norm(U3hatnew-U3hatold, "F")>eps)  & iter<maxiter)
  {iter=iter+1
  U1hatold<-U1hatnew
  U2hatold<-U2hatnew
  U3hatold<-U3hatnew
  iter=iter+1
  tt2=U3hatold%x%U1hatold
  for (k in 1:K){U2int[,,k]=XU2[,,k]%*%solve(tt2)%*%aperm(XU2[,,k], perm=c(2,1))}
  U2hatnew<-apply(U2int, MARGIN=c(1,2), sum)/(n3*n1*K)
  tt3=U2hatnew%x%U1hatold
  for (k in 1:K){U3int[,,k]=XU3[,,k]%*%solve(tt3)%*%aperm(XU3[,,k], perm=c(2,1))}
  U3hatnew<-apply(U3int, MARGIN=c(1,2), sum)/(n2*n1*K)
  tt1=U3hatnew%x%U2hatnew
  for (k in 1:K){U1int[,,k]=XU1[,,k]%*%solve(tt1)%*%aperm(XU1[,,k], perm=c(2,1))}
  U1hatnew<-apply(U1int, MARGIN=c(1,2), sum)/(n2*n3*K)
  }
  if(iter==maxiter) {
    Iter<-c("Did not converge at maximum number of iterations given eps. Try to increase maxit and/or decrease eps.")
    Convergence<-FALSE
  } else  {
    Iter<-iter
    Convergence<-TRUE}
  Shat<-U3hatnew%x%U2hatnew%x%U1hatnew
  #Printing out results
  list(Call=formula,
       Convergence=Convergence,
       Iter=Iter,
       Xmeanhat=Xmean,
       First=c(names(mf)[2], levels(mf[,2])),
       U1hat=U1hatnew,
       Standardized.U1hat=U1hatnew/U1hatnew[1,1],
       Second=c(names(mf)[3], levels(mf[,3])),
       U2hat=U2hatnew,
       Standardized.U2hat=U2hatnew*(U1hatnew[1,1]),
       Third=c(names(mf)[4], levels(mf[,4])),
       U3hat=U3hatnew,
       Shat=Shat
    )
}

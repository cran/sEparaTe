% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mle2d_svc.R
\name{mle2d_svc}
\alias{mle2d_svc}
\title{Maximum likelihood estimation of the parameters of a matrix normal distribution}
\usage{
mle2d_svc(value2d, Id1, Id2, subject, data_2d, eps, maxiter, startmat)
}
\arguments{
\item{value2d}{from the formula value2d ~ Id1 + Id2}

\item{Id1}{from the formula value2d ~ Id1 + Id2}

\item{Id2}{from the formula value2d ~ Id1 + Id2}

\item{subject}{the replicate, also called individual}

\item{data_2d}{the name of the matrix data}

\item{eps}{the threshold in the stopping criterion for the iterative mle algorithm}

\item{maxiter}{the maximum number of iterations for the iterative mle algorithm}

\item{startmat}{the value of the second factor variance-covariance matrix used for
initialization, i.e., to start the algorithm and obtain the initial estimate
of the first factor variance-covariance matrix}
}
\description{
Maximum likelihood estimation for the parameters of a matrix normal
distribution \bold{X}, which is characterized by a simply separable
variance-covariance structure. In the general case, which is the case
considered here, two unstructured factor
variance-covariance matrices determine the covariability of random
matrix entries, depending on the row (one factor matrix)
and the column (the other factor matrix) where two
\bold{X}-entries are. In the required function, the Id1 and Id2 variables correspond to
the row and column subscripts, respectively; \dQuote{value2d}
indicates the observed variable.
}
\section{Output}{


\dQuote{Convergence}, TRUE or FALSE

\dQuote{Iter}, will indicate the number of iterations needed for the mle algorithm to converge

\dQuote{Xmeanhat}, the estimated mean matrix (i.e., the sample mean)

\dQuote{First}, the row subscript, or the second column in the data file

\dQuote{U1hat}, the estimated variance-covariance matrix for the rows

\dQuote{Standardized.U1hat}, the standardized estimated variance-covariance matrix
for the rows; the standardization is performed by dividing each entry of U1hat
by entry(1, 1) of U1hat

\dQuote{Second}, the column subscript, or the third column in the data file

\dQuote{U2hat}, the estimated variance-covariance matrix for the columns

\dQuote{Standardized.U2hat}, the standardized estimated variance-covariance
matrix for the columns; the standardization is performed by multiplying each
entry of U2hat by entry(1, 1) of U1hat

\dQuote{Shat}, is the sample variance-covariance matrix computed from of the vectorized data matrices
}

\section{References}{


Dutilleul P. 1990. Apport en analyse spectrale d'un periodogramme
modifie et modelisation des series chronologiques avec repetitions en vue de leur comparaison
en frequence. D.Sc. Dissertation, Universite catholique
de Louvain, Departement de mathematique.

Dutilleul P. 1999. The mle algorithm for the matrix
normal distribution. Journal of Statistical Computation
and Simulation 64: 105-123.
}

\examples{
output <- mle2d_svc(data2d$value2d, data2d$Id1, data2d$Id2, data2d$K, data_2d = data2d)
output

}

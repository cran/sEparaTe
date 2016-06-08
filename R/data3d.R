#' Three dimensional data set
#'
#' An i.i.d. random sample of size 11 from a 2 x 2 x 2
#'    tensor normal distribution, for a small numerical example of the use of the
#'    functions mle3d_svc and lrt3d_svc from the \pkg{sEparaTe} package
#'
#' @format A frame (excluding the headings) with 88 lines of data and 5 variables:
#' \describe{
#'   \item{K}{an integer ranging from 1 to 11, the size of an i.i.d. random sample from a 2 x 2 x 2 tensor matrix normal distribution}
#'   \item{Id3}{an integer ranging from 1 to 2, the number of rows of the 3rd-order tensor normal distribution}
#'   \item{Id4}{an integer ranging from 1 to 2, the number of columns of the 3rd-order tensor normal distribution}
#'   \item{Id5}{an integer ranging from 1 to 2, the number of edges of the 3rd-order tensor normal distribution}
#'   \item{value3d}{the sample data for the observed variable}
#' }

"data3d"

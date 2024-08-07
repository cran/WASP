% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn.R
\name{knn}
\alias{knn}
\title{Modified k-nearest neighbour conditional bootstrap or regression function estimation with extrapolation}
\usage{
knn(
  x,
  z,
  zout,
  k = 0,
  pw,
  reg = TRUE,
  nensemble = 100,
  tailcorrection = TRUE,
  tailprob = 0.25,
  tailfac = 0.2,
  extrap = TRUE
)
}
\arguments{
\item{x}{A vector of response.}

\item{z}{A matrix of existing predictors.}

\item{zout}{A matrix of predictor values the response is to be estimated at.}

\item{k}{The number of nearest neighbours used. The default value is 0, indicating Lall and Sharma default is used.}

\item{pw}{A vector of partial weights of the same length of z.}

\item{reg}{A logical operator to inform whether a conditional expectation should be output or not nensemble, Used if reg=F and represents the number of realisations that are generated Value.}

\item{nensemble}{An integer the specifies the number of ensembles used. The default is 100.}

\item{tailcorrection}{A logical value, T (default) or F, that denotes whether a reduced value of k (number of nearest neighbours) should be used in the tails of any conditioning plane. Whether one is in the tails or not is determined based on the nearest neighbour response value.}

\item{tailprob}{A scalar that denotes the p-value of the cdf (on either extreme) the tailcorrection takes effect. The default value is 0.25.}

\item{tailfac}{A scalar that specifies the lowest fraction of the default k that can be used in the tails. Depending on the how extreme one is in the tails, the actual k decreases linearly from k (for a p-value greater than tailprob) to tailfac*k proportional to the actual p-value of the nearest neighbour response, divided by tailprob. The default value is 0.2.}

\item{extrap}{A logical value, T (default) or F, that denotes whether a kernel extraplation method is used to predict x.}
}
\value{
A matrix of responses having same rows as zout if reg=T, or having nensemble columns is reg=F.
}
\description{
Modified k-nearest neighbour conditional bootstrap or regression function estimation with extrapolation
}
\examples{
\donttest{
# AR9 model   x(i)=0.3*x(i-1)-0.6*x(i-4)-0.5*x(i-9)+eps
data.ar9 <- data.gen.ar9(500)
x <- data.ar9$x # response
z <- data.ar9$dp # possible predictors

zout <- ts(data.gen.ar9(500, ndim = ncol(z))$dp) # new input

xhat1 <- xhat2 <- x
xhat1 <- knn(x, z, zout, k = 5, reg = TRUE, extrap = FALSE) # without extrapolation
xhat2 <- knn(x, z, zout, k = 5, reg = TRUE, extrap = TRUE) # with extrapolation

ts.plot(ts(x), ts(xhat1), ts(xhat2), col = c("black", "red", "blue"),
ylim = c(-5, 5), lwd = c(2, 2, 1))
}
}
\references{
Sharma, A., Tarboton, D.G. and Lall, U., 1997. Streamflow simulation: A nonparametric approach. Water resources research, 33(2), pp.291-308.

Sharma, A. and O'Neill, R., 2002. A nonparametric approach for representing interannual dependence in monthly streamflow sequences. Water resources research, 38(7), pp.5-1.
}

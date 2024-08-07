% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modwt_vt.R
\name{modwt.vt}
\alias{modwt.vt}
\title{Variance Transformation Operation - MODWT}
\usage{
modwt.vt(
  data,
  wf,
  J,
  boundary,
  cov.opt = "auto",
  flag = "biased",
  detrend = FALSE,
  backward = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{data}{A list of response x and dependent variables dp.}

\item{wf}{Name of the wavelet filter to use in the decomposition.}

\item{J}{Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).}

\item{boundary}{Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.}

\item{cov.opt}{Options of Covariance matrix sign. Use "pos", "neg", or "auto".}

\item{flag}{Biased or Unbiased variance transformation, c("biased","unbiased").}

\item{detrend}{Detrend the input time series or just center, default (F).}

\item{backward}{Detrend the input time series or just center, default (F).}

\item{verbose}{A logical indicating if some “progress report” should be given.}
}
\value{
A list of 8 elements: wf, J, boundary, x (data), dp (data), dp.n (variance transformed dp), and S (covariance matrix).
}
\description{
Variance Transformation Operation - MODWT
}
\examples{
### real-world example
data(Ind_AWAP.2.5)
data(obs.mon)
data(SPI.12)
x <- window(SPI.12, start = c(1950, 1), end = c(2009, 12))
dp <- window(obs.mon, start = c(1950, 1), end = c(2009, 12))

op <- par(mfrow = c(ncol(dp), 1), pty = "m", mar = c(1, 4, 1, 2))
for (id in sample(Ind_AWAP.2.5, 1)) {
  data <- list(x = x[, id], dp = dp)
  dwt <- modwt.vt(data, wf = "d4", J = 7, boundary = "periodic", cov.opt = "auto")

  for (i in 1:ncol(dp)) {
    ts.plot(dwt$dp[, i], dwt$dp.n[, i], xlab = NA, col = c("black", "red"), lwd = c(2, 1))
  }
}
par(op)

### synthetic example
# frequency, sampled from a given range
fd <- c(3, 5, 10, 15, 25, 30, 55, 70, 95)

data.SW1 <- data.gen.SW(nobs = 512, fp = 25, fd = fd)
dwt.SW1 <- modwt.vt(data.SW1, wf = "d4", J = 7, boundary = "periodic", cov.opt = "auto")

x.modwt <- waveslim::modwt(dwt.SW1$x, wf = "d4", n.levels = 7, boundary = "periodic")
dp.modwt <- waveslim::modwt(dwt.SW1$dp[, 1], wf = "d4", n.levels = 7, boundary = "periodic")
dp.vt.modwt <- waveslim::modwt(dwt.SW1$dp.n[, 1], wf = "d4", n.levels = 7, boundary = "periodic")

sum(sapply(dp.modwt, var))
var(dwt.SW1$dp[, 1])
sum(sapply(dp.vt.modwt, var))
var(dwt.SW1$dp.n[, 1])

data <- rbind(
  sapply(dp.modwt, var) / sum(sapply(dp.modwt, var)),
  sapply(dp.vt.modwt, var) / sum(sapply(dp.vt.modwt, var))
)

bar <- barplot(data, beside = TRUE, col = c("red", "blue"))
lines(x = bar[2, ], y = sapply(x.modwt, var) / sum(sapply(x.modwt, var)))
points(x = bar[2, ], y = sapply(x.modwt, var) / sum(sapply(x.modwt, var)))

}
\references{
Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962.

Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2020). A wavelet-based tool to modulate variance in predictors: an application to predicting drought anomalies. Environmental Modelling & Software, 135, 104907.
}

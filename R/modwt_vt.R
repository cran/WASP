#--------------------------------------------------------------------------------
#' Variance Transformation Operation - MODWT
#' @param data		  A list of response x and dependent variables dp.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt   Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#' @param flag      Biased or Unbiased variance transformation, c("biased","unbiased").
#' @param detrend   Detrend the input time series or just center, default (F).
#' @param backward  Detrend the input time series or just center, default (F).
#' @param verbose	  A logical indicating if some “progress report” should be given.
#'
#' @return A list of 8 elements: wf, J, boundary, x (data), dp (data), dp.n (variance transformed dp), and S (covariance matrix).
#' @import waveslim
#' @export
#'
#' @references Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962.
#' @references Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2020). A wavelet-based tool to modulate variance in predictors: an application to predicting drought anomalies. Environmental Modelling & Software, 135, 104907.
#'
#' @examples
#' ### real-world example
#' data(Ind_AWAP.2.5)
#' data(obs.mon)
#' data(SPI.12)
#' x <- window(SPI.12, start = c(1950, 1), end = c(2009, 12))
#' dp <- window(obs.mon, start = c(1950, 1), end = c(2009, 12))
#'
#' op <- par(mfrow = c(ncol(dp), 1), pty = "m", mar = c(1, 4, 1, 2))
#' for (id in sample(Ind_AWAP.2.5, 1)) {
#'   data <- list(x = x[, id], dp = dp)
#'   dwt <- modwt.vt(data, wf = "d4", J = 7, boundary = "periodic", cov.opt = "auto")
#'
#'   for (i in 1:ncol(dp)) {
#'     ts.plot(dwt$dp[, i], dwt$dp.n[, i], xlab = NA, col = c("black", "red"), lwd = c(2, 1))
#'   }
#' }
#' par(op)
#'
#' ### synthetic example
#' # frequency, sampled from a given range
#' fd <- c(3, 5, 10, 15, 25, 30, 55, 70, 95)
#'
#' data.SW1 <- data.gen.SW(nobs = 512, fp = 25, fd = fd)
#' dwt.SW1 <- modwt.vt(data.SW1, wf = "d4", J = 7, boundary = "periodic", cov.opt = "auto")
#'
#' x.modwt <- waveslim::modwt(dwt.SW1$x, wf = "d4", n.levels = 7, boundary = "periodic")
#' dp.modwt <- waveslim::modwt(dwt.SW1$dp[, 1], wf = "d4", n.levels = 7, boundary = "periodic")
#' dp.vt.modwt <- waveslim::modwt(dwt.SW1$dp.n[, 1], wf = "d4", n.levels = 7, boundary = "periodic")
#'
#' sum(sapply(dp.modwt, var))
#' var(dwt.SW1$dp[, 1])
#' sum(sapply(dp.vt.modwt, var))
#' var(dwt.SW1$dp.n[, 1])
#'
#' data <- rbind(
#'   sapply(dp.modwt, var) / sum(sapply(dp.modwt, var)),
#'   sapply(dp.vt.modwt, var) / sum(sapply(dp.vt.modwt, var))
#' )
#'
#' bar <- barplot(data, beside = TRUE, col = c("red", "blue"))
#' lines(x = bar[2, ], y = sapply(x.modwt, var) / sum(sapply(x.modwt, var)))
#' points(x = bar[2, ], y = sapply(x.modwt, var) / sum(sapply(x.modwt, var)))
#'
modwt.vt <- function(data, wf, J, boundary, cov.opt = "auto",
                     flag = "biased", detrend = FALSE, backward=FALSE, verbose=TRUE) {
  # initialization
  x <- data$x
  dp <- as.matrix(data$dp)
  mu.dp <- apply(dp, 2, mean)
  
  # reverse data to normal timeline
  if(backward) {
    x <- rev(x)
    dp <- apply(dp, 2, rev)
  }
  
  # output 
  ndim <- ncol(dp)
  n <- nrow(dp)
  S <- matrix(nrow = J + 1, ncol = ndim)
  dp.n <- matrix(nrow = n, ncol = ndim)
  modwt.dp <- vector("list", ndim)
  Wn.list <- vector('list', NCOL(dp))

  for (i in 1:ndim) {
    # center or detrend
    if (!detrend) {
      dp.c <- scale(dp[, i], scale = F)
    } else {
      # dp.c <- lm(dp[,i]~c(1:n))$resid
      dp.c <- dp[, i] - smooth.spline(1:n, dp[, i], spar = 1)$y
    }

    # MODWT - variance decomposition
    modwt.dp[[i]] <- waveslim::modwt(dp.c, wf = wf, n.levels = J, boundary = boundary)
    B <- matrix(unlist(modwt.dp[[i]]), ncol = J + 1, byrow = FALSE)

    Bn <- scale(B)[1:n, ]
    Wn.list[[i]] <- Bn
    V <- as.numeric(apply(B, 2, sd))

	if(verbose){
		# dif <- sum(abs(imodwt(modwt.dp[[i]]) - dp.c)) # this is equivalent to MODWT-MRA
		dif <- sum(abs(Bn %*% V - dp.c))
		if (dif > 10^-10) print(paste0("Difference between reconstructed and original series: ", dif))
	}

    # variance transformation
    cov <- cov(x, Bn[seq_along(x), ])
    # cat("Biased: ", round(cov,3),"\n")

    if (flag == "unbiased") { ### unbiased wavelet variance - only change cov
      modwt.dp.n <- non.bdy(modwt.dp[[i]], wf = wf, method = "modwt")

      B.n <- matrix(unlist(modwt.dp.n), ncol = J + 1, byrow = FALSE)
      cov <- cov(x, scale(B.n)[seq_along(x), ], use = "pairwise.complete.obs")

      # if unbiased cov is not available then use biased
      cov[is.na(cov)] <- cov(x, Bn[seq_along(x), ])[is.na(cov)]
      # cat("Unbiased: ",round(cov,3),"\n")
    }

    if (cov.opt == "pos") cov <- cov else if (cov.opt == "neg") cov <- -cov
    S[, i] <- as.vector(cov)

    Vr <- as.numeric(cov / norm(cov, type = "2") * sd(dp.c))
    if (!detrend) {
      dp.n[, i] <- Bn %*% Vr + mu.dp[i]
    } else {
      # dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
      dp.n[, i] <- Bn %*% Vr + smooth.spline(1:n, dp[, i], spar = 1)$y
    }

    # check the correlation after vt then decide the direction of C
    if (cov.opt == "auto") {
      if (cor.test(dp.n[, i], dp[, i])$estimate < 0 & cor.test(dp.n[, i], dp[, i])$p.value < 0.05) cov <- -cov
      S[, i] <- as.vector(cov)

      Vr <- as.numeric(cov / norm(cov, type = "2") * sd(dp.c))
      if (!detrend) {
        dp.n[, i] <- Bn %*% Vr + mu.dp[i]
      } else {
        # dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
        dp.n[, i] <- Bn %*% Vr + smooth.spline(1:n, dp[, i], spar = 1)$y
      }
    }

    # dif.var <- abs(var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    # if(dif.var>0.15) print(paste0("Variance difference between transformed
    #                  and original series by percentage: ", dif.var*100))
  }

  # reverse data to normal timeline
  if(backward) {
    x <- rev(x)
    dp <- apply(dp, 2, rev)
    dp.n <- apply(dp.n, 2, rev)
  }
  
  dwt <- list(
    wavelet = wf,
    J = J,
    boundary = boundary,
    x = x,
    dp = dp,
    dp.n = dp.n,
    S = S,
    Wn = Wn.list
  )
  class(dwt) <- "modwt"

  return(dwt)
}

#--------------------------------------------------------------------------------
#' Variance Transformation Operation for Validation
#'
#' @param data		  A list of response x and dependent variables dp.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt       A class of "modwt" data. Output from modwt.vt().
#' @param detrend   Detrend the input time series or just center, default (F).
#' @param backward  Detrend the input time series or just center, default (F).
#' @param verbose	  A logical indicating if some “progress report” should be given.
#'
#' @return A list of 8 elements: wf, J, boundary, x (data), dp (data), dp.n (variance transformed dp), and S (covariance matrix).
#' @export
#'
#' @references Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962. doi:10.1029/2019wr026962
#'
#' @examples
#' data(rain.mon)
#' data(obs.mon)
#'
#' ## response SPI - calibration
#' SPI.cal <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(1979,12)),sc=12)
#' #SPI.cal <- SPEI::spi(window(rain.mon, start = c(1949, 1), end = c(1979, 12)), scale = 12)$fitted
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for (id in 1:ncol(SPI.cal)) {
#'   x <- window(SPI.cal[, id], start = c(1950, 1), end = c(1979, 12))
#'   dp <- window(obs.mon, start = c(1950, 1), end = c(1979, 12))
#'   data.list[[id]] <- list(x = as.numeric(x), dp = matrix(dp, nrow = nrow(dp)))
#' }
#'
#' ## variance transformation - calibration
#' dwt.list <- lapply(data.list, function(x) {
#'   modwt.vt(x, wf = "d4", J = 7, boundary = "periodic", cov.opt = "auto")
#' })
#'
#' ## response SPI - validation
#' SPI.val <- SPI.calc(window(rain.mon, start=c(1979,1), end=c(2009,12)),sc=12)
#' #SPI.val <- SPEI::spi(window(rain.mon, start = c(1979, 1), end = c(2009, 12)), scale = 12)$fitted
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for (id in 1:ncol(SPI.val)) {
#'   x <- window(SPI.val[, id], start = c(1980, 1), end = c(2009, 12))
#'   dp <- window(obs.mon, start = c(1980, 1), end = c(2009, 12))
#'   data.list[[id]] <- list(x = as.numeric(x), dp = matrix(dp, nrow = nrow(dp)))
#' }
#'
#' # variance transformation - validation
#' dwt.list.val <- lapply(
#'   seq_along(data.list),
#'   function(i) modwt.vt.val(data.list[[i]], J = 7, dwt.list[[i]])
#' )
#'
#' ## plot original and reconstrcuted predictors for each station
#' for (i in seq_along(dwt.list.val)) {
#'   # extract data
#'   dwt <- dwt.list.val[[i]]
#'   x <- dwt$x # response
#'   dp <- dwt$dp # original predictors
#'   dp.n <- dwt$dp.n # variance transformed predictors
#'
#'   plot.ts(cbind(x, dp))
#'   plot.ts(cbind(x, dp.n))
#' }
modwt.vt.val <- function(data, J, dwt, detrend = FALSE, backward=FALSE, verbose=TRUE) {

  # initialization
  x <- data$x
  dp <- as.matrix(data$dp)
  wf <- dwt$wavelet
  boundary <- dwt$boundary
  mu.dp <- apply(dp, 2, mean)
  
  # reverse data
  if(backward) {
    x <- rev(x)
    dp <- apply(dp, 2, rev)
  }
  
  # output
  ndim <- ncol(dp)
  n <- nrow(dp)
  dp.n <- matrix(nrow = n, ncol = ndim)
  modwt.dp <- vector("list", ndim)
  for (i in 1:ndim) {
    # center or detrend
    if (!detrend) {
      dp.c <- scale(dp[, i], scale = F)
    } else {
      # dp.c <- lm(dp[,i]~c(1:n))$resid
      dp.c <- dp[, i] - smooth.spline(1:n, dp[, i], spar = 1)$y
    }

    # MODWT - variance decomposition
    modwt.dp[[i]] <- waveslim::modwt(dp.c, wf = wf, n.levels = J, boundary = boundary)
    B <- matrix(unlist(modwt.dp[[i]]), ncol = J + 1, byrow = FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B, 2, sd))

	if(verbose){
		# dif <- sum(abs(imodwt(modwt.dp[[i]]) - dp.c)) # this is equivalent to MODWT-MRA
		dif <- sum(abs(Bn %*% V - dp.c))
		if (dif > 10^-10) print(paste0("Difference between reconstructed and original series: ", dif))
	}

    # in case different J
    cov <- rep(0, J + 1)
    if (length(dwt$S[, i]) > (J + 1)) {
      cov <- dwt$S[, i][1:(J + 1)]
    } else {
      cov[seq_along(dwt$S[, i])] <- dwt$S[, i]
    }

    Vr <- as.numeric(cov / norm(cov, type = "2") * sd(dp.c))

    if (!detrend) {
      dp.n[, i] <- Bn %*% Vr + mu.dp[i]
    } else {
      # dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
      dp.n[, i] <- Bn %*% Vr + smooth.spline(1:n, dp[, i], spar = 1)$y
    }

    # dif.var <- abs(var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    # if (dif.var > 0.15) print(paste0("Variance difference between transformed
    #                     and original series by percentage: ", dif.var * 100))
  }

  # reverse data to normal timeline
  if(backward) {
    x <- rev(x)
    dp <- apply(dp, 2, rev)
    dp.n <- apply(dp.n, 2, rev)
  }
  
  dwt <- list(
    wavelet = wf,
    J = J,
    boundary = boundary,
    x = x,
    dp = dp,
    dp.n = dp.n,
    S = dwt$S
  )
  class(dwt) <- "modwt"

  return(dwt)
}

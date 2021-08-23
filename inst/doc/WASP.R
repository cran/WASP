## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)
library(rmarkdown)
library(knitr)
library(kableExtra)

## ----setup--------------------------------------------------------------------
library(WASP)
library(ggplot2)
library(fracdiff)

library(FNN)
library(synthesis)
library(waveslim)
library(cowplot)

## ----wavelet-transforms-------------------------------------------------------
  # data generation
  x <- arima.sim(list(order = c(1,0,0), ar = 0.6), n = 512)
  #x <- as.numeric(scale(data.gen.Rossler(time = seq(0, 50, length.out = 512))$x, scale=F))

  #for(wf in c("mb4","w4","bs3.1")){ #not working
  for(wf in c("haar","d4","d8","d16", "fk4","la8","bl14")){ #working ones
    print(paste0("Wavelet filter: ", wf))
    #----------------------------------------------------------------------------
    # wavelet family, extension mode and package
    #wf <- "haar" # wavelet family D8 or db4
    boundary <- "periodic"
    if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1

    #Maximum decomposition level J
    n <- length(x)
    J <- ceiling(log(n/(2*v-1))/log(2))-1 #(Kaiser, 1994)
    
    cov = rnorm(J+1, sd=2); Vr <- as.numeric(cov/norm(cov,type="2")*sd(x))
    #----------------------------------------------------------------------------
    #DWT-MRA
    print("-----------DWT-MRA-----------")
    x.mra <- waveslim::mra(x,wf=wf, J=J, method="dwt", boundary=boundary)
    x.mra.m <- matrix(unlist(x.mra), ncol=J+1)
    
    x.n <- scale(x.mra.m)%*%Vr; var(x.n)-var(x)

    message(paste0("Additive decompostion: ", isTRUE(all.equal(as.numeric(x),rowSums(x.mra.m)))))
    message(paste0("Variance decompostion: ", isTRUE(all.equal(var(x),sum(apply(x.mra.m,2,var))))))

    #----------------------------------------------------------------------------
    #MODWT
    print("-----------MODWT-----------")
    x.modwt <- waveslim::modwt(x, wf=wf, n.levels = J, boundary = boundary)
    x.modwt.m <- matrix(unlist(x.modwt), ncol=J+1)
    
    x.n <- scale(x.modwt.m)%*%Vr; var(x.n)-var(x)
    
    message(paste0("Additive decompostion: ", isTRUE(all.equal(as.numeric(x),rowSums(x.modwt.m)))))
    message(paste0("Variance decompostion: ", isTRUE(all.equal(var(x),sum(apply(x.modwt.m,2,var))))))

    #----------------------------------------------------------------------------
    #a trous
    print("-----------AT-----------")
    x.at <- at.wd(x, wf=wf, J=J, boundary=boundary)
    x.at.m <- matrix(unlist(x.at), ncol=J+1)
    
    # x.mra.modwt <- waveslim::mra(x,wf=wf, J=J, method="modwt", boundary=boundary)
    # x.mra.modwt <- matrix(unlist(x.mra.modwt), ncol=J+1)
    # 
    # print(sum(abs(x.at.m-x.mra.modwt)))

    message(paste0("Additive decompostion: ", isTRUE(all.equal(as.numeric(x),rowSums(x.at.m)))))
    message(paste0("Variance decompostion: ", isTRUE(all.equal(var(x),sum(apply(x.at.m,2,var))))))

    if(isTRUE(all.equal(x.at.m,x.modwt.m))) message(paste0("AT and MODWT is equivalent using ", wf))

  }
  

## ----wt-overview--------------------------------------------------------------
# data generation
x <- arima.sim(list(order = c(1,0,0), ar = 0.6), n = 512)

if(TRUE){
  #----------------------------------------------------------------------------
  # wavelet family, extension mode and package
  wf <- "la8" # wavelet family D8 or db4
  boundary <- "periodic"
  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
  
  #Maximum decomposition level J
  n <- length(x)
  J <- ceiling(log(n/(2*v-1))/log(2))-1 #(Kaiser, 1994)
  
  limits.x<- c(0,n); limits.y <- NULL #c(-3,3)
  #----------------------------------------------------------------------------
  #DWT-MRA
  x.mra <- waveslim::mra(x,wf=wf, J=J, method="dwt", boundary=boundary)
  x.mra.m <- matrix(unlist(x.mra), ncol=J+1)
  
  if(TRUE){
    mra.plot(x, x.mra.m, limits.x, limits.y, ylab="X", col="red", type="details", main=paste0("DWT","(",wf,")"))
    p1 <- recordPlot()
  }
  
  #----------------------------------------------------------------------------
  #MODWT
  x.modwt <- waveslim::modwt(x, wf=wf, n.levels = J, boundary = boundary)
  x.modwt.m <- matrix(unlist(x.modwt), ncol=J+1)
  
  if(TRUE){
    mra.plot(x, x.modwt.m, limits.x, limits.y, ylab="X", col="red", type="coefs", main=paste0("MODWT","(",wf,")"))
    p2 <- recordPlot()
  }
  
  #----------------------------------------------------------------------------
  #a trous
  x.at <- at.wd(x, wf=wf, J=J, boundary=boundary)
  x.at.m <- matrix(unlist(x.at), ncol=J+1)
  
  if(TRUE){
    mra.plot(x, x.at.m, limits.x, limits.y, ylab="X", col="red", type="coefs", main=paste0("AT","(",wf,")"))
    p3 <- recordPlot()
  }

  #----------------------------------------------------------------------------
  #plot and save
  cowplot::plot_grid(p1,p2,p3, ncol=3, labels = c("(a)","(b)","(c)"), label_size = 12)
}

## ----wt-alignment-------------------------------------------------------------
# data generation
nobs=64
set.seed(2020)
#x <- data.gen.ar4(nobs)$dp[,1]
x <- WASP::data.gen.SW(nobs,fp=5, fd=5)$x
#plot.ts(x)

if(TRUE){
  #----------------------------------------------------------------------------
  # wavelet family, extension mode and package
  wf <- "haar" # wavelet family D8 or db4
  boundary <- "periodic"
  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
  
  #Maximum decomposition level J
  n <- length(x)
  J <- ceiling(log(n/(2*v-1))/log(2))-1 #(Kaiser, 1994)
  J <- 3
  
  #----------------------------------------------------------------------------
  #DWT-MRA
  x.mra <- waveslim::mra(x,wf=wf, J=J, method="dwt", boundary=boundary)
  x.mra.m <- matrix(unlist(x.mra), ncol=J+1)
  
  #----------------------------------------------------------------------------
  #MODWT
  x.modwt <- waveslim::modwt(x, wf=wf, n.levels = J, boundary = boundary)
  x.modwt.m <- matrix(unlist(x.modwt), ncol=J+1)
  
  #----------------------------------------------------------------------------
  #a trous
  x.at <- at.wd(x, wf=wf, J=J, boundary=boundary)
  x.at.m <- matrix(unlist(x.at), ncol=J+1)
  
  ylab = c("D", "A")
  limits.x<- c(0,n); limits.y <- c(-4,4); cols <- c("black","red","blue")
  if(TRUE){
    op <- par(mfcol = c(J, 1), mar = c(0, 3, 0.2, 1), oma = c(2, 
        1, 1, 1), mgp = c(1, 1, 0), bg = "transparent", cex.lab=1.5,
        pty = "m", ps = 8)
    for (i in 1:J) {
        ts.plot(as.ts(cbind(x, x.mra.m[,i],x.modwt.m[,i])), gpars = list(axes = FALSE, xlab = NA, ylab = paste0(ylab[1], i), xaxs = "i", xlim = limits.x, ylim = limits.y, col=cols))
        box()
        legend("topleft", legend=c("Orignal", "DWT-MAR", "MODWT/AT"), col=c("black","red","blue"), lty=1, cex=1.45)
    }
    axis(side = 1, at = seq(limits.x[1], limits.x[2], by = 4), 
        labels = seq(limits.x[1], limits.x[2], by = 4),cex.axis=2)
    legend("topleft", legend=c("Orignal", "DWT-MAR", "MODWT/AT"), col=c("black","red","blue"), lty=1, cex=1.45)
    par(op)
  }
}

## ----wt-decomposition-level---------------------------------------------------
sample <- seq(10, by=100, length.out=20)
v=2
tmp <- NULL
for(n in sample){
J1 <- floor(log(n/(2*v-1))/log(2));J #(Kaiser, 1994)
J2 <- floor(log2(n/(2*v-1)-1));J #Cornish, C. R., Bretherton, C. S., & Percival, D. B. (2006)
J3 <- floor(log10(n));J #(Nourani et al., 2008)

tmp <- cbind(tmp, c(J1,J2,J3))
}
print(tmp)

## ----wavelet-variance---------------------------------------------------------
if(TRUE){
  # Generates simulated long-memory time series data from the fractional ARIMA(p,d,q) model.
  nobs = 360 
  set.seed(2020)
  x <- fracdiff.sim(nobs, d = .4)$series

  #x <- rnorm(nobs)
  
  #x <- arima.sim(list(order = c(0,1,0)), n = nobs-1)
  
  # f=c(16,32,64,128); t <- seq(0,1,length.out = nobs); const=c(0.01,0.1,1,10)
  # x.m <- sapply(1:length(f), function(i) const[i]*sin(2*pi*f[i]*t)); apply(x.m,2,var)
  # x <- rowSums(x.m) 
  
  #x <- rnorm(nobs)
  
  #x <- arima.sim(list(order = c(2,0,0), ar=c(0.75,-0.5)), n = nobs)
  
  # plot.ts(x)
  # acf(x)
} else {
  data("rain.mon"); data("obs.mon")
  spi = window(SPEI::spi(rain.mon[,5], scale=12)$fitted, start=c(1950,1), end=c(2009,12))
  x <- as.numeric(spi); plot.ts(x)
  acf(x)  
  
  x <- obs.mon[,1]
  nobs=length(x)
}

wf="haar"; n=nobs
if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
J <- floor(log(n/(2*v-1))/log(2))
boundary = "periodic"
xaxis <- 1/2^((J-1):0)
#----------------------------------------------------------
###dwt-mra
x.mra <- waveslim::mra(padding(x,"zero"), wf=wf, J=J, method="dwt", boundary=boundary)
x.mra.m <- matrix(unlist(x.mra), ncol=J+1)
v.mra <- apply(x.mra.m[1:nobs,],2, var)

# analysis of variance
x.mra.bw <- non.bdy(x.mra, wf, method="mra") #replacing boundary wavelet coefficients with NA.
x.mra.bw.m <- matrix(unlist(lapply(x.mra.bw,function(i) i[1:n])), ncol=J+1)
x.dwt.var <- wave.variance(lapply(x.mra.bw,function(i) i[1:n]), type="eta3")

Bn <- scale(x.mra.bw.m)

cov <- cov(x, Bn[1:length(x),], use="pairwise.complete.obs");cov
Vr <- as.numeric(cov/norm(cov,type="2")*sd(x))

Bn[is.na(Bn)]=0
x1 <- Bn%*%Vr
x2 <- scale(x.mra.m[1:nobs,])%*%Vr

sum(abs(x1-x))

ts.plot(cbind(x,x1,x2), col=1:3)

#----------------------------------------------------------
###modwt
x.modwt <- waveslim::modwt(x, wf=wf, n.levels = J, boundary = boundary)
x.modwt.m <- matrix(unlist(x.modwt), ncol=J+1)
v.modwt <- apply(x.modwt.m, 2, var)

# analysis of variance
x.modwt.bw <- brick.wall(x.modwt, wf, method="modwt")
x.modwt.bw.m <- matrix(unlist(x.modwt.bw), ncol=J+1)
x.modwt.var <- wave.variance(x.modwt.bw, type="eta3")

#----------------------------------------------------------
###plot
ylim = c(10^floor(log10(min(v.mra,v.modwt,na.omit(x.modwt.var[,1]),na.omit(x.dwt.var[,1])))), ceiling(max(v.mra,v.modwt,na.omit(x.modwt.var[,1]),na.omit(x.dwt.var[,1]))))
op <- par(mfrow=c(2,1), las=1, mar=c(5,4,2,2)+.1)
if(TRUE){
matplot(xaxis, x.dwt.var[-(J+1),1], type="b", log="xy",
        xaxt="n", pch=1, lty=1, col=c(1,4), ylim=ylim, 
        xlab="Wavelet Scale", ylab="", lwd=2)
matlines(xaxis, x.dwt.var[-(J+1),2:3], lty=2, col="grey")
matlines(xaxis, v.mra[-(J+1)], type="b",
         pch=2, lty=1, col=4)
axis(side=1, at=xaxis)
legend("topright", c("Wavelet variance", "Actual"),
       lty=1, col=c(1,4), bty="n")

matplot(xaxis, x.modwt.var[-(J+1),1], type="b", log="xy",
        xaxt="n", pch=1, lty=1, col=c(1,4), ylim=ylim, 
        xlab="Wavelet Scale", ylab="", lwd=2)
matlines(xaxis, x.modwt.var[-(J+1),2:3], lty=2, col="grey")
matlines(xaxis, v.modwt[-(J+1)], type="b",
         pch=2, lty=1, col=4)
axis(side=1, at=xaxis)
legend("topright", c("Wavelet variance", "Actual"),
       lty=1, col=c(1,4), bty="n")

}
par(op)

## ----variance-transform, warning=TRUE-----------------------------------------
#-------------------------------------------------------------------
op <- par()
wf <- "d16"
flag = switch(1, "biased", "unbiased")
sample=10000
ylim <- c(-55,55)
if(TRUE){
  set.seed(2020)
  ###synthetic example - Rossler

  s=0.1
  ts.list <- list()
  for(i in seq_along(s)){
    ts.r <- data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2), time = seq(0, 50, length.out = sample))
  
    #add noise
    ts.r$x <- ts(ts.r$x + rnorm(n = sample, mean=0, sd=s[i]))
    ts.r$y <- ts(ts.r$y + rnorm(n = sample, mean=0, sd=s[i]))
    ts.r$z <- ts(ts.r$z + rnorm(n = sample, mean=0, sd=s[i]))
  
    ts.list[[i]]<- ts.r
  }
  
  data.list <- lapply(ts.list, function(ts) list(x=ts$z, dp=cbind(ts$x,ts$y)))
  
  lab.names <- c("x","y")

} else {
  
  ###Real-world example
  data("obs.mon"); data("rain.mon")


  SPI.12 <- SPEI::spi(rain.mon[,5],scale=12)$fitted
  x <- window(SPI.12,start=c(1950,1),end=c(2009,12))
  dp <- window(obs.mon,start=c(1950,1),end=c(2009,12))
  
  data.list <- list(list(x=x,dp=dp)) 
  sample=length(x)
  
  lab.names <- colnames(obs.mon)
  ylim <- NULL
  
}

#-------------------------------------------------------------------
p.list <- list(); dp.list <- list()
if(wf!="haar") mode.opts <- c("MRA", "MODWT","AT")[1:3]
if(wf=="haar") mode.opts <- c("MRA", "MODWT","AT")[1:2]
for(mode in mode.opts){

  cov.opt <- switch(1,"auto","pos","neg")
  if(mode=="MRA") method <- switch(1,"dwt","modwt")

  # wavelet family, extension mode and package
  #wf <- switch(mode, "MRA"="haar", "MODWT"="haar", "AT"="haar")
  pad="zero"
  boundary <- "periodic"
  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1

  #Maximum decomposition level J
  n <- sample
  J <- ceiling(log(n/(2*v-1))/log(2)) - 1 #(Kaiser, 1994)
  #J <- floor(log(n/(2*v-1))/log(2))

  #variance transfrom - calibration
  if(mode=="MRA"){
    dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt, flag))
  } else if(mode=="MODWT") {
    dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, boundary, cov.opt, flag))
  } else {
    dwt.list<- lapply(data.list, function(x) at.vt(x, wf, J, boundary, cov.opt, flag))
  }

  for(j in 1:length(dwt.list)){
    dwt <- dwt.list[[j]]

    par(mfrow=c(ncol(dwt$dp),1),mar=c(0,2.5,2,1),
        oma = c(2, 1, 0, 0), # move plot to the right and up
        mgp=c(1.5, 0.5, 0), # move axis labels closer to axis
        pty="m",bg = "transparent",
        ps=6)

    #plot(dwt$x, type="l", xlab=NA, ylab="SPI12", col="red")
    #plot(dwt$x, type="l", xlab=NA, ylab="Rain", col="red")
    for(i in 1:ncol(dwt$dp))
      ts.plot(cbind(dwt$dp[,i],dwt$dp.n[,i]),xlab=NA,ylab=paste0(lab.names[i]),ylim=ylim,
              col=c("black","blue"),lwd=c(1,2))

    p.list[[length(p.list)+1]] <- recordPlot()

    dp.list[[length(dp.list)+1]] <- dwt$dp.n

  }
}

# if(!isTRUE(all.equal(dp.list[[1]],dp.list[[3]]))) warning(paste0("DWT-MRA and MODWT is not equivalent using ", wf))
# if(isTRUE(all.equal(dp.list[[2]],dp.list[[3]]))) warning(paste0("AT and MODWT is equivalent using ", wf))
#----------------------------------------------------------------------------
#plot and save
fig <- cowplot::plot_grid(plotlist =p.list, nrow=1, labels= c("(a)","(b)","(c)"))
fig

par(op)

## ----haar---------------------------------------------------------------------
if(TRUE){
  # Generates simulated long-memory time series data from the fractional ARIMA(p,d,q) model.
  nobs = 1024
  x <- fracdiff.sim(nobs, d = .4)$series
  #acf(x)
  #fd=25
  #x <- WASP::data.gen.SW(nobs, fd=25, sd.y=0)$x
  #plot.ts(x)
  
  x <- data.gen.Rossler(time=seq(0,50, length.out = nobs))$z

} else {
  data("rain.mon"); data("obs.mon")
  spi = window(SPEI::spi(rain.mon[,5], scale=12)$fitted, start=c(1950,1), end=c(2009,12))
  x <- as.numeric(spi); plot.ts(x)
  #acf(x)  
  i=7
  x <- obs.mon[,i];print(colnames(obs.mon)[i])
  
  data("slf.ts")
  x <- slf.ts
  
  nobs=length(x)
}


J=7; boundary = "periodic"
xaxis <- 1/2^((J-1):0)

for(wf in c("haar","d4","la8")){
  print(paste0("----------------",wf,"-----------------------"))
  #----------------------------------------------------------
  ###dwt-mra
  x.mra <- waveslim::mra(padding(x,"zero"), wf=wf, J=J, method="dwt", boundary=boundary)
  x.mra.m <- matrix(unlist(x.mra), ncol=J+1)
  
  # overview
  plot.ts(x.mra.m[1:nobs,])
  
  print(round(cor(x.mra.m[1:nobs,]),3))
  
  
  #----------------------------------------------------------
  ###modwt
  x.modwt <- waveslim::modwt(x, wf=wf, n.levels = J, boundary = boundary)
  x.modwt.m <- matrix(unlist(x.modwt), ncol=J+1)
  
  # overview
  plot.ts(x.modwt.m)
  
  print(round(cor(x.modwt.m),3))
  
  # #----------------------------------------------------------
  # ###alignment effects
  # par(mfrow=c((J+1),1), mar=rep(1,4))
  # for(i in 1:(J+1)){
  # 
  #   sub <- which(as.integer(time(x))==1986)
  #   ts.plot(cbind(x[sub],x.modwt.m[sub, i], x.mra.m[sub,i]), col=c("black","red","blue"))
  # 
  # }

}


## ----wavelet-transforms-compare-----------------------------------------------
  # data generation
  #x <- arima.sim(list(order = c(1,0,0), ar = 0.6), n = 512)
  data <- data.gen.ar1(nobs=512)
  x=data$dp[,data$true.cpy[1]]; y=data$x

  plot.ts(cbind(x,y))
  RMSE <- NULL
  #for(wf in c("mb4","w4","bs3.1")){ #not working
  for(wf in c("haar","d4","d8","d16", "fk4","la8","bl14")){ #working ones
    print(paste0("Wavelet filter: ", wf))
    #----------------------------------------------------------------------------
    # wavelet family, extension mode and package
    #wf <- "haar" # wavelet family D8 or db4
    boundary <- "periodic"
    if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1

    #Maximum decomposition level J
    n <- length(x)
    if(wf=="haar") J<-ceiling(log(n/(2*v-1))/log(2))-1 else J<-ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
    
    cov = rnorm(J+1, sd=2); Vr <- as.numeric(cov/norm(cov,type="2")*sd(x))
    #----------------------------------------------------------------------------
    #DWT-MRA
    print("-----------DWT-MRA-----------")
    x.mra <- waveslim::mra(x,wf=wf, J=J, method="dwt", boundary=boundary)
    x.mra.m <- matrix(unlist(x.mra), ncol=J+1)
    
    m1 <- lm(y~x.mra.m);
    rmse1 <- sqrt(mean(m1$residuals^2))
    
    #----------------------------------------------------------------------------
    #MODWT
    print("-----------MODWT-----------")
    x.modwt <- waveslim::modwt(x, wf=wf, n.levels = J, boundary = boundary)
    x.modwt.m <- matrix(unlist(x.modwt), ncol=J+1)
    
    x.n <- scale(x.modwt.m)%*%Vr; var(x.n)-var(x)
    
    m2 <- lm(y~x.modwt.m)
    rmse2 <- sqrt(mean(m2$residuals^2))

    #----------------------------------------------------------------------------
    #a trous
    print("-----------AT-----------")
    x.at <- at.wd(x, wf=wf, J=J, boundary=boundary)
    x.at.m <- matrix(unlist(x.at), ncol=J+1)
    
    x.n <- scale(x.at.m)%*%Vr; var(x.n)-var(x)

    m3 <- lm(y~x.at.m)
    rmse3 <- sqrt(mean(m3$residuals^2))
    
    RMSE <- rbind(RMSE, c(rmse1,rmse2,rmse3))

  }
  
print(RMSE)
  

## ----optimal-variance-transformation, warning=TRUE----------------------------
if(TRUE){
###Synthetic example
#data generation
set.seed(2020)
sample = 512
#frequency, sampled from a given range
fd <- c(3,5,10,15,25,30,55,70,95)
#data <- WASP::data.gen.SW(nobs=sample,fp=25,fd=fd)
data <- WASP::data.gen.SW(nobs=sample,fp=c(15,25,30),fd=fd)

# ts = data.gen.Rossler(time = seq(0, 50, length.out = sample))
# data <- list(x=ts$z, dp=cbind(ts$x, ts$y))

} else {
###Real-world example
data("obs.mon"); data("rain.mon")

if(1){ #SPI12 as response
  SPI.12 <- SPEI::spi(rain.mon[,5],scale=12)$fitted
  x <- window(SPI.12,start=c(1950,1),end=c(2009,12))
  dp <- window(obs.mon,start=c(1950,1),end=c(2009,12))
  
} else {#rainfall as response
  x <- window(rain.mon[,5],start=c(1950,1),end=c(2009,12))
  dp <- window(obs.mon,start=c(1950,1),end=c(2009,12))
}
data <- list(x=x,dp=dp)
sample=length(x)

}

#plot.ts(cbind(data$x,data$dp))

tab.list <- list()
mode.opts <- c("MRA", "MODWT","AT")
for(mode in mode.opts){
  print(mode)
  
  #cov.opt <- switch(2,"auto","pos","neg")
  if(mode=="MRA") method <- switch(1,"dwt","modwt")
  
  # wavelet family, extension mode and package
  #wf <- switch(mode, "MRA"="haar", "MODWT"="haar", "AT"="haar")
  wf="haar"
  pad="zero"
  boundary <- "periodic"
  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
  
  #Maximum decomposition level J
  n <- sample
  J <- ceiling(log(n/(2*v-1))/log(2))-1 #(Kaiser, 1994)
  
  tab <- NULL
  for(cov.opt in c("auto","pos","neg")){
    # variance transfrom - calibration
    if(mode=="MRA"){
      dwt<- dwt.vt(data, wf, J, method, pad, boundary, cov.opt)
    } else if(mode=="MODWT") {
      dwt<- modwt.vt(data, wf, J, boundary, cov.opt)
    } else {
      dwt<- at.vt(data, wf, J, boundary, cov.opt)
    }
    
    # optimal prediction accuracy
    opti.rmse <- NULL
    dp.RMSE <- NULL; dp.n.RMSE <- NULL
    S <- dwt$S; ndim=ncol(S)
    for(i in 1:ndim){
      x <- dwt$x
      dp <- dwt$dp[,i]
      dp.n <- dwt$dp.n[,i]
      
      #ts.plot(cbind(dp,dp.n), col=1:2)
  
      dp.RMSE <- c(dp.RMSE, sqrt(mean(lm(x~dp)$residuals^2)))
      dp.n.RMSE <- c(dp.n.RMSE, sqrt(mean(lm(x~dp.n)$residuals^2)))
  
      #small difference due to the reconstruction
      opti.rmse <- c(opti.rmse, sqrt((n-1)/n*(var(x)-sum(S[,i]^2)*var(dp)/var(dp.n)))) 
      #opti.rmse <- c(opti.rmse, sqrt((n-1)/n*(var(x)-sum(S[,i]^2))))
    }
  
  tab <- rbind(tab, cbind(dp.RMSE, dp.n.RMSE, opti.rmse))
  }
  
  rownames(tab) <- rep(c("auto","pos","neg"),each=ndim)
  tab.list[[length(tab.list)+1]] <- tab
} 

print(tab.list)


## ----high-order-transformation------------------------------------------------
#-------------------------------------------------------------------
###Real-world example
data("obs.mon")
data("rain.mon")
op <- par()
station.id = 5
lab.names <- colnames(obs.mon)[c(1,3,4,5,7)]

if(TRUE){ #SPI12 as response
  SPI.12 <- SPEI::spi(rain.mon,scale=12)$fitted
  x <- window(SPI.12,start=c(1950,1),end=c(2009,12))
  dp <- window(obs.mon[,lab.names],start=c(1950,1),end=c(2009,12))
  
} else {#rainfall as response
  x <- window(rain.mon,start=c(1950,1),end=c(2009,12))
  dp <- window(obs.mon[,lab.names],start=c(1950,1),end=c(2009,12))
  
}

data.list <- lapply(station.id, function(id) list(x=x[,id],dp=dp))


ylim=data.frame(GPH=c(700,900),TDP700=c(5,25),TDP500=c(5,25),EPT=c(300,330),
                UWND=c(-5,25),VWND=c(-5,10),MSLP=c(-1, 1))[c(1,3,4,5,7)]

#-------------------------------------------------------------------
p.list <- list()
RMSE <- NULL
mode.opts <- c("MRA", "MODWT","AT")[1:2]
for(mode in mode.opts){
  
  cov.opt <- switch(1,"auto","pos","neg")
  if(mode=="MRA") method <- switch(1,"dwt","modwt")
  
  # wavelet family, extension mode and package
  wf <- switch(mode, "MRA"="d4", "MODWT"="haar", "AT"="haar")
  pad="zero"
  boundary <- "periodic"
  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
  
  # Maximum decomposition level J
  n <- nrow(x)
  J <- ceiling(log(n/(2*v-1))/log(2)) - 1 #(Kaiser, 1994)
  
  # high order variance transfromation
  dwt.list <- lapply(data.list, function(data) stepwise.VT(data, mode=mode, wf=wf))
  
  for(j in seq_len(length(dwt.list))){
    dwt <- dwt.list[[j]]
    cpy <- dwt$cpy
    
    MSE <- NULL
    for(i in seq_len(length(cpy))){
      m1 <- sqrt(FNN::knn.reg(train=dwt$dp[,1:i], y=dwt$x)$PRESS/n)
      m2 <- sqrt(FNN::knn.reg(train=dwt$dp.n[,1:i], y=dwt$x)$PRESS/n)
      
      MSE <- rbind(MSE, c(m1,m2))
    }
    
    RMSE <- cbind(RMSE, MSE)
    
    par(mfrow=c(length(cpy),1),mar=c(0,4,2,1),
        oma = c(2, 1, 0, 0), # move plot to the right and up
        mgp=c(1.5, 0.5, 0), # move axis labels closer to axis
        pty="m",bg = "transparent",
        ps=8)
    
    #plot(dwt$x, type="l", xlab=NA, ylab="SPI12", ylim=c(-3,3),col="red")
    #plot(dwt$x, type="l", xlab=NA, ylab="Rain", col="red")
    for(i in seq_len(length(cpy))){
      ts.plot(dwt$dp[,i],dwt$dp.n[,i],xlab=NA,ylab=paste0(lab.names[cpy[i]]), #ylim=ylim[,i],
              col=c("black","blue"),lwd=c(1,2))
    }
    
    p.list[[length(p.list)+1]] <- recordPlot()
    
  }
}

#-------------------------------------------------------------------
#plot and save
cowplot::plot_grid(plotlist =p.list, nrow=1, labels= c("(a)","(b)","(c)"))
par(op)
#-------------------------------------------------------------------
#RMSE when more predictors are included
tab1 <- round(RMSE,3)
tab1 <- cbind(1:nrow(tab1), tab1)
colnames(tab1) <- c("No. of Predictors", rep(c("Original","Transformed"), length(mode.opts)))

kable(tab1, caption = "Comparison of prediction accuracy using both original and transformed high order predictors", booktabs = T) %>%
  kable_styling(latex_options = c("HOLD_position"),position = "center") %>%
#  add_header_above(c(" " = 1, "DWT-MRA" = 2, "MODWT" = 2, "AT" = 2))
  add_header_above(c(" " = 1, "DWT-MRA" = 2, "MODWT/AT" = 2)) 


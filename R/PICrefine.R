gaussfit = function(x,y){
  library(stats)
  a1 <- max(y[2:(length(y)-1)])
  b1 <- x[which(y==max(y[2:(length(y)-1)]))]
  fit <- nls(y~a1*exp(-((x-b1)/c1)^2),start=list(a1=a1,b1=b1,c1=1),
             control=list(maxiter=500,tol=1e-3,minFactor=1/512,printEval=F,warnOnly=T))
  yp <- predict(fit)
  sse <- sum((yp-y)^2)
  sst <- sum((y-mean(y))^2)
  r.square <- 1-sse/sst
  if (r.square<0){r.square <- 0}
  return(r.square)
}

get.sharpness = function(sig){
  p <- which(sig[2:(length(sig)-1)]==max(sig[2:(length(sig)-1)]))+1
  sharpness <- 0
  for (i in 2:p){
    if (sig[i-1]>0.1){temp <- ((sig[i]-sig[i-1]))/sig[i-1]
    }else{temp <- 0}
    sharpness <- sharpness+temp
  }
  for (j in p:length(sig)-1){
    if (sig[i+1]>0.1){temp <- ((sig[i]-sig[i+1])/sig[i+1])
    }else{temp <- 0}
    sharpness <- sharpness+temp
  }
  return(sharpness)
}

waveft = function(omega,scales)
{
  StpFrq <- omega[2]
  NbFrq  <- length(omega)
  SqrtNbFrq <- sqrt(NbFrq)
  cfsNORM <- sqrt(StpFrq)*SqrtNbFrq
  NbSc <- length(scales)
  wft <- matrix(0,NbSc,NbFrq)
  
  mul <- sqrt(scales/gamma(2+0.5))*cfsNORM
  for (jj in 1:NbSc){
    scapowered = (scales[jj]*omega)
    expnt = -(scapowered^2)/2
    wft[jj,] = mul[jj]*(scapowered^2)*exp(expnt)
  }
  return (wft)
}

cwtft = function(sig)
{
  library(stats)
  options(warn =-1)
  val <- sig
  meanSIG <- mean(val)
  xx <- val-meanSIG
  n <- length(xx)
  dt <- 1
  
  omega <- (1:floor(n/2))
  omega <- omega*((2*pi)/(n*dt))
  omega <- c(0, omega,-(omega[seq(floor((n-1)/2),1,-1)]))
  f <- fft(xx)
  
  # set scales
  s0 <- 2*dt
  ds <- 0.4875
  NbSc <- floor(log2(n)/ds)
  scales <- s0*2^(0:(NbSc-1)*ds)
  param <- 2
  
  psift <- waveft(omega,scales)
  mat <- matrix(0,NbSc,length(f))
  for (ii in 1:NbSc){
    mat[ii,] <- f
  }
  cwtcfs <- t(mvfft(t(mat*psift),inverse = TRUE)/ncol(mat))
  cwtcfs <- cwtcfs[,1:n]
  result <- matrix(0,nrow(cwtcfs),ncol(cwtcfs))
  for (jj in 1:nrow(cwtcfs)){
    result[jj,] <- as.numeric(cwtcfs[jj,])
  }
  result <- t(result)
  colnames(result) <- scales
  return (result)
}

subfun = function(x,PICs,tInl)
{
  PIC <- PICs$PICs[[x]]
  rtmin <- as.numeric(min(PIC[,1]))
  rtmax <- as.numeric(max(PIC[,1]))
  inte <- PIC[,2]
  
  rt1 <- seq(from=rtmin,to=rtmax,by=tInl)
  inte1 <- approx(x=PIC[,1],y=PIC[,2],rt1)$y
  wCoefs <- cwtft(inte1)
  noise <- as.numeric(quantile(abs(wCoefs[,1]),0.8))
  siglevel <- max(wCoefs)
  
  r.square <- try(gaussfit(rt1,inte1),T)
  if(is.numeric(r.square)==F){r.square<-0.0}
  sharpness <- get.sharpness(inte1)
  SNR <- siglevel/noise
  
  peakInfo <- list()
  peakInfo$r.square <- 1.5*r.square
  peakInfo$sharpness <- sharpness
  peakInfo$SNR <- SNR
  
  return(peakInfo)
}

PICrefine = function(PICs,n=1)
{
  library(parallel)
  library(doParallel)
  library(foreach)
  subfun <- subfun
  cwtft <- cwtft

  output <- list()
  tInl <- PICs$rt[2]-PICs$rt[1]

  cl <- makeCluster(getOption("cl.cores", detectCores(logical = F)))
  registerDoParallel(cl)
  peakInfo <- foreach(x=1:length(PICs$PICs)) %dopar%
    subfun(x,PICs,tInl)
  stopCluster(cl)
  scores <- matrix(0,length(peakInfo),4)
  for(jj in 1:length(peakInfo)){
    scores[jj,1] <- peakInfo[[jj]]$r.square
    scores[jj,2] <- peakInfo[[jj]]$sharpness
    if(peakInfo[[jj]]$SNR>3){
      scores[jj,3]<-1}else{scores[jj,3]<-peakInfo[[jj]]$SNR/3}
  }
  scores[,2] <- ((scores[,2]-mean(scores[,2]))/var(scores[,2])+1)*0.4
  scores[,4] <- rowSums(scores[,1:3])
  colnames(scores) <- c("Gaussian similarity","sharpness","SNR","total")
  refine <- which(scores[,4]<mean(scores[,4])-n*var(scores[,4]))
  PIClist <- PICs$PICs
  PIClist[refine] <-NULL
  ind <- setdiff(1:nrow(scores),refine)
  pInfo <- cbind(PICs$Info,scores[,4])[ind,]
  colnames(pInfo) <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","Intensity","score")

  output$PICs <- PIClist
  output$Info <- pInfo
  output$scores <- scores[ind,]
  output$rt <- PICs$rt
  return (output)
}
smoothPIC <- function(PIC,lambda){
  PICs <- PIC$PICs
  corr <- matrix(0,length(PICs),2)
  colnames(corr) <- c("corr1","corr2")
  for (i in 1:length(PICs)){
    sig <- PICs[[i]]
    # freq <- PIC$DataInfo$freq
    # rt <- seq(sig[1,1],sig[nrow(sig),1],freq)
    rt0 <- PIC$rt
    rt <- rt0[rt0>=sig[1,1]&rt0<=sig[nrow(sig),1]]
    inte <- approx(sig[,1],sig[,2],xout=rt)$y
    sig <- cbind(rt,inte)
    colnames(sig) <- c('rt','intensity')
    try({result <- gaussfit(sig[,1],sig[,2],lambda,rt0)
    corr[i,] <- result$correlation
    PICs[[i]] <- result$sig}, silent=T)
  }
  Info <- cbind(PIC$Info,corr)
  PIC$smoothedPIC <- PICs
  PIC$Info <-Info
  return(PIC)
}


smoother <- function(y,lambda){
  M <- length(y)
  E <- diag(M)
  D <- diff(E)
  z <- solve((E+lambda*t(D)%*%D),y)
  return(z)
}

gaussfit = function(x,y,lambda,rt0){
  library(stats)
  a1 <- max(y[2:(length(y)-1)])
  b1 <- x[which(y==max(y[2:(length(y)-1)]))]
  fit <- nls(y~a1*exp(-((x-b1)/c1)^2),start=list(a1=a1,b1=b1,c1=1),
             control=list(maxiter=500,tol=1e-3,minFactor=1/512,printEval=F))
  coeff <- summary(fit)$coefficients
  a1 <- coeff['a1',1]
  b1 <- coeff['b1',1]
  c1 <- coeff['c1',1]
  sigma <- c1/sqrt(2)
  freq <- x[2]-x[1]
  num <- min(round(1.5*sigma/freq),round(length(x)/3))
  x1 <- rt0[rt0>=(min(x)-num*freq)&rt0<=(max(x)+num*freq)]
  y1 <- a1*exp(-((x1-b1)/c1)^2)
  y <- c(rep(0,num),y,rep(0,num))
  corr1 <- cor(y,y1)
  y <- smoother(y,lambda)
  corr2 <- cor(y,y1)
  correlation <- c(corr1,corr2)
  return(list(sig=cbind(x1,y),correlation=correlation))
}

getSpectrum <- function(X,R,mztol=0.05,pw=20){
  rt1 <- X$rt
  rt2 <- R$rt
  rt1 <- round(rt1,2)
  rt2 <- round(rt2,2)
  if(length(rt1)!=length(rt2)){
    rt1 <- rt1[rt1>ceiling(max(rt1[1],rt2[1]))&rt1<floor(min(rt1[length(rt1)],rt2[length(rt2)]))]
    rt2 <- rt2[rt2>ceiling(max(rt1[1],rt2[1]))&rt2<floor(min(rt1[length(rt1)],rt2[length(rt2)]))]
  }
  if(length(rt1)!=length(rt2)){
    rt1 <- rt1[1:min(length(rt1),length(rt2))]
    rt2 <- rt2[1:min(length(rt1),length(rt2))]
  }
  x <- rep(0,length(rt1))
  r <- rep(0,length(rt2))
  w <- 0
  X.id <- which(X$Info[,'corr1']>0.95&X$Info[,'corr2']>0.98)
  R.id <- which(R$Info[,'corr1']>0.95&R$Info[,'corr2']>0.98)
  index <- seq(1,length(x),round(pw/(rt1[2]-rt1[1])))
  for (i in 1:floor(length(index)/2)){
    peak1 <- X.id[X$Info[X.id,'rt']>rt1[index[2*i-1]]&X$Info[X.id,'rt']<rt1[index[2*i]]]
    peak2 <- R.id #[R$Info[R.id,'rt']>rt1[index[2*i-1]]&R$Info[R.id,'rt']<rt1[index[2*i]]]
    r.peak <- NULL
    if (length(peak1)>0){
      temp <- order(X$Info[peak1,'Intensity'])
      for (j in 1:length(peak1)){
        x.peak <- peak1[which(temp==j)]
        r.candidate <- peak2[abs(R$Info[peak2,'mz']-X$Info[x.peak,'mz'])<mztol]
        if (length(r.candidate)==1){
          r.peak <- r.candidate
          break
        }
      }
    }
    if (length(r.peak)!=1){next}
    x.PIC <- X$smoothedPIC[[x.peak]]
    r.PIC <- R$smoothedPIC[[r.peak]]
    w <- max(w,nrow(x.PIC),nrow(r.PIC))

    x[which(rt1%in%round(x.PIC[,1],2))] <- x[which(rt1%in%round(x.PIC[,1],2))]+x.PIC[,2]
    r[which(rt2%in%round(r.PIC[,1],2))] <- r[which(rt2%in%round(r.PIC[,1],2))]+r.PIC[,2]
  }
  w <- round(1.2*w)
  return(list(x=x,r=r,w=w,x.rt=rt1,r.rt=rt2))
}

mwfft <- function(x,r,w){
  p1 <- t(as.matrix(x))
  p2 <- t(as.matrix(r))
  shifts <- rep(0,length(p1))
  p1e <- cbind(t(rep(0,w)),p1,t(rep(0,w)))
  p2e <- cbind(t(rep(0,w)),p2,t(rep(0,w)))
  
  shifts[1] <- FFTcorr(t(as.matrix(p1e[1,1:w])), t(as.matrix(p2e[1,1:w])),w)
  mat <- rbind(matrix(0,w-1,1),shifts[1])
  
  for (i in 2:(length(p1e)-w)){
    inds <- seq(i,(i+w-1),1)
    shifts[i] <- FFTcorr(t(as.matrix(p1e[1,inds])), t(as.matrix(p2e[1,inds])),w)
    z <- c(mat[2:w,ncol(mat)],shifts[i])
    mat <- cbind(mat,z)
  }
  colnames(mat) <- NULL
  shifts <- shifts[w+1:length(shifts)]
  mat <- mat[,(w+1):ncol(mat)]
  shift_mode <- rep(0,ncol(mat))
  for (ii in 1:length(shift_mode)){
    shift_mode[ii] <- as.numeric(names(table(mat[,ii]))[which.max(table(mat[,ii]))])
  }
  shifts_refine <- shift_mode
  dlr <- diff(shifts_refine)
  inds <- which(dlr!=0)
  
  if (length(inds)>0){
    for (i in 2:(length(inds))){
      s <- shifts_refine[inds[i]+1]-shifts_refine[inds[i]]
      if (i!=length(inds)&&((2*abs(s)+1)>=(inds[i]-inds[i-1])||(2*abs(s)+1)>=(inds[i+1]-inds[i]))){
        shifts_refine[inds[i]:inds[i+1]] <- shifts_refine[inds[i]]
      }else{if(i==length(inds)&&(2*abs(s)+1)>=(inds[i]-inds[i-1])){
        shifts_refine[inds[i]:length(shifts_refine)]= shifts_refine[inds[i]]
      }}
    }
  }
  
  dlr <- diff(shifts_refine)
  inds <- which(dlr!=0)
  if (length(inds)>0){
    s <- shifts_refine[inds[1]+1]-shifts_refine[inds[1]]
    p1_aligned <- as.numeric(p1[1:(inds[1]-abs(s)-1)])
    if (shifts_refine[1]>0){
      p1_aligned <- c(matrix(0,1,shifts_refine[1]),p1_aligned)
    }else{if(shifts_refine[1]<0){
      p1_aligned <- p1_aligned[abs(shifts_refine[1])+1:length(p1_aligned)]
    }}
    
    for (i in 1:length(inds)){
      s <- shifts_refine[inds[i]+1]-shifts_refine[inds[i]]
      if ((inds[i]-abs(s))<=0){
        ss=1:(inds[i]+abs(s))
      }else{if((inds[i]+abs(s))>length(shifts_refine)){
        ss=((inds[i]-abs(s)):length(shifts_refine))
      }else{ss=((inds[i]-abs(s)):(inds[i]+abs(s)))
      }}
      sv <- p1[ss]
      si <- order(sv)
      v <- sv[si]
      sv <- as.numeric(sv)
      if (s<0){
        sv <- sv[-(si[1:abs(s)])]
      }else{if(length(sv)>=si[1]+1){sv <- c(sv[1:si[1]],as.numeric(v[1]*rep(1,abs(s))),sv[si[1]+1:length(sv)])
      }else{sv <- c(sv[1:si[1]],as.numeric(v[1]*rep(1,abs(s))))}}
      sv <- sv[is.na(sv)==F]
      p1_aligned <- c(p1_aligned,sv)
      if (i<length(inds)){
        sn <- shifts_refine[inds[i+1]+1]-shifts_refine[inds[i+1]]
        sns <- 0
        try(sns <- seq((inds[i]+abs(s)+1),(inds[i+1]-abs(sn)-1),1),T)
        snv <- as.numeric(p1[sns])
      }else{
        sn <- 0-shifts_refine[inds[i]+1]
        sns <- 0
        try(sns <- seq(inds[i]+abs(s)+1,length(p1),1),T)
        snv <- as.numeric(p1[sns])
        if (length(p1_aligned)>length(x)){p1_aligned <- p1_aligned(1:length(x))
        }else{if(sn>0){snv <- c(snv,snv[length(snv)]*matrix(1,1,sn))}else{
          if(sn<0){snv <- snv[-(snv[length(snv)-abs(sn)+1:length(snv)])]}}}
      }
      snv <- as.numeric(snv)
      p1_aligned <- c(p1_aligned,snv)
    }
  }else{p1_aligned <- p1}
  x_aligned <- p1_aligned[1:length(x)]
  return(list(x_aligned=x_aligned,shifts_refine=shifts_refine))
}

FFTcorr <- function(spectrum, target, shift){
  M <- ncol(target)
  diff <- 1000000
  curdiff <- 2^(1:20)
  for (i in 1:20){
    curdiff <- ((2^i)-M)
    if (curdiff>0&curdiff<diff){diff <- curdiff}
  }
  target <- cbind(target,t(rep(0,diff)))
  spectrum <- cbind(spectrum,t(rep(0,diff)))
  M <- M+diff
  X <- fft(as.numeric(target))
  Y <- fft(as.numeric(spectrum))
  R <- X*Conj(Y)
  R <- R/M
  rev <- fft(R,inverse=T)/length(rev)
  vals <- Re(rev)
  maxpos <- 1
  maxi <- -1
  if (M<shift){shift <- M}
  for (i in 1:shift){
    if (vals[i] > maxi){
      maxi = vals[i]
      maxpos = i
    }
    if (vals[length(vals)-i+1] > maxi){
      maxi = vals[length(vals)-i+1];
      maxpos = length(vals)-i+1;
    }
  }
  if (maxi < 0.1){lag <- 0
  return(lag)}
  if (maxpos > length(vals)/2){
    lag = maxpos-length(vals)-1
    return(lag)
  }else{lag <- maxpos-1
  return(lag)}
}

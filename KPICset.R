setClass("KPICSet", representation(peakmat="matrix",rt="list",sample="vector",PICset="list"))
         
KPICset <- function(files,range=10,level=500,alpha=0.2,Iterate=4000,n=3,refine=NA){
  library(parallel)
  library(iterators)
  library(foreach)
  library(doParallel)
  library(Ckmeans.1d.dp)
  output <- new("KPICSet")
  filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                   "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
  filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  info <- file.info(files)
  listed <- list.files(files[info$isdir], pattern = filepattern,
                       recursive = TRUE, full.names = TRUE)
  sample <- list.files(files[info$isdir], pattern = filepattern,
                       recursive = F, full.names = F)
  for (i in 1:length(sample)){
    temp <- unlist(strsplit(sample[i],'[.]'))
    sample[[i]] <- temp[1]
  }
  files <- c(files[!info$isdir], listed)
  result <- getPIC(files[1],range=range,level=level,alpha=alpha,Iterate=Iterate,n=n)
  if (is.na(refine)==F){
    result <- PICrefine(result,refine)
  }
  peakmat <- matrix(1,nrow(result$Info),1)
  colnames(peakmat) <- "sample"
  peakmat <- cbind(result$Info,peakmat)
  output@PICset$PICset[[1]] <- result
  output@rt[[1]] <- result$rt
  cat(1,'/',length(files), "is done","\n")
  for (i in 2:length(files)){
    result <- getPIC(files[i],range=range,level=level,alpha=alpha,Iterate=Iterate,n=n)
    if (is.na(refine)==F){
      result <- PICrefine(result,refine)
    }
    peakmat1 <- matrix(i,nrow(result$Info),1)
    colnames(peakmat1) <- "sample"
    peakmat1 <- cbind(result$Info,peakmat1)
    peakmat <- rbind(peakmat,peakmat1)
    output@PICset$PICset[[i]] <- result
    output@rt[[i]] <- result$rt
    cat(i,'/',length(files), "is done","\n")
  }
  output@peakmat <- peakmat
  output@sample <- sample
  return(output)
}

rtcor <- function(xset,mztol=0.05,pw=20,lambda=1.5){
  R <- xset@PICset$PICset[[1]]
  R <- smoothPIC(R,lambda)
  peaks <- xset@peakmat
  rt_cor <- round(peaks[,'rt'],2)
  peaks <- cbind(peaks,rt_cor)
  for (nsam in 2:length(xset@sample)){
    X <- xset@PICset$PICset[[nsam]]
    X <- smoothPIC(X,lambda)
    Spectrum <- getSpectrum(X,R,mztol=mztol,pw=pw)
    align.result <- mwfft(Spectrum$x,Spectrum$r,Spectrum$w)
    x.rt <- x.align <- round(X$rt,2)
    x.align[x.rt%in%Spectrum$x.rt] <- x.rt[which(x.rt%in%Spectrum$x.rt)+align.result$shifts_refine]
    rt.raw <- X$Info[,'rt']
    rt.align <- matrix(rt.raw,length(rt.raw),1)
    colnames(rt.align) <- 'rt'
    for (i in 1:length(rt.raw)){
      rt.align[i,] <- x.align[which(abs(x.rt-rt.raw[i])==min(abs(x.rt-rt.raw[i])))]
    }
    X$Info <- cbind(X$Info,rt.align)
    peaks[which(peaks[,'sample']==nsam),'rt_cor'] <- rt.align
    xset@PICset$PICset[[nsam]] <- X
    cat(nsam,'/',length(xset@sample), "is done","\n")
  }
  xset@PICset$PICset[[1]] <- R
  xset@peakmat <- peaks
  return(xset)
}
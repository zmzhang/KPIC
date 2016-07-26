setClass("KPICSet", representation(peakmat="matrix",rt="list",sample="vector",n_features="numeric",PICset="list"))

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
  n_features <- rep(0,length(sample))
  peakmat <- matrix(1,nrow(result$Info),1)
  colnames(peakmat) <- "sample"
  peakmat <- cbind(result$Info,peakmat)
  output@PICset <- list(result)
  output@rt[[1]] <- result$rt
  n_features[1] <- nrow(peakmat)
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
    output@PICset <- c(output@PICset,list(result))
    output@rt[[i]] <- result$rt
    n_features[i] <- nrow(peakmat1)
    cat(i,'/',length(files), "is done","\n")
  }
  output@peakmat <- peakmat
  output@sample <- sample
  output@n_features <- n_features
  return(output)
}

rtcor <- function(xset,files,mzSegSize=.5,shift=10,lambda=1.5,ref=NA){
  lambda <- max(0,lambda)
  if (is.na(ref)){
    TICs <- getTIC(files,method='BPC')
    ref <- TICs$ref
  }
  peakmat <- xset@peakmat
  shift <- round(shift/(xset@rt[[1]][2]-xset@rt[[1]][1]))
  mzlist1 <- seq((min(xset@peakmat[,'mz'])-0.5*mzSegSize),
                (max(xset@peakmat[,'mz'])+0.5*mzSegSize),mzSegSize)
  mzlist2 <- mzlist1+0.5*mzSegSize
  peakmat[,'rt'] <- rt_cor <- round(peakmat[,'rt'],2)
  peakmat <- cbind(peakmat,rt_cor)
  xset@peakmat <- peakmat
  for (i in 1:(length(mzlist1)-1)){
    xset <- seg.rtcor(xset,ref,mzlist1[i],mzlist1[i+1],shift=shift,lambda=lambda)
    xset <- seg.rtcor(xset,ref,mzlist2[i],mzlist2[i+1],shift=shift,lambda=lambda)
    if (i%%20==0){cat(round(i/(length(mzlist1)-1)*100),'%',' is done','\n')}
  }
  return(xset)
}
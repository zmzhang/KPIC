LoadData = function(fileName)
{
  library (mzR)
  rawMs <- openMSfile(fileName)
  MsPeak <- peaks(rawMs)
  Info <- runInfo(rawMs)
  freq <- (Info$dEndTime-Info$dStartTime)/(Info$scanCount-1)
  Info$freq <- freq
  
  num <- sapply(MsPeak, length)/2
  Mat <- cbind(rep(1:length(MsPeak), num),do.call(rbind, MsPeak))
  Mat[,1] <- (Mat[,1]-1)*freq+Info$dStartTime
  
  data <- list()
  data$Mat <-Mat
  data$Info <- Info
  return (data)
}

subPIC = function(mat,range,level,alpha=0.3,Iterate=4000,tInl,gap)
{
  library(Ckmeans.1d.dp)
  PICs <- list()
  Info <- matrix(0,0,ncol=7)
  
  mat <- mat[order(mat[,3],decreasing=T),]
  id <- seq(1,nrow(mat),1)
  options(warn =-1)
  nPIC <- 1
  for (i in 1:Iterate){
    if(length(id)==0){break}
    rInte <- mat[id[1],3]
    if (rInte<level){
      break
    }
    rTime <- mat[id[1],1]
    rMz <- mat[id[1],2]
    ssid <- id[which(abs(mat[id,1]-rTime)<range)]
    mzdiff <- abs(mat[ssid,2]-rMz)
    jj <- mzdiff<alpha*range/tInl*gap
    ssid <- ssid[jj]
    mzdiff <- mzdiff[jj]
    output <- Ckmeans.1d.dp(mzdiff^2, k=c(1,5))
    mincenter <- min(output$centers)
    tClu <- which(output$centers==mincenter)
    ssid <- ssid[which(output$cluster==tClu)]
    
    if(length(ssid)>5){
      PImat <- cbind(ssid,mat[ssid,])
      PImat <- PImat[order(PImat[,2]),]
      pTime <- PImat[,2]
      if(length(pTime)!=length(unique(pTime))){
        pTime <- unique(PImat[,2])
        for (j in 1:length(pTime)){
          aa <- which(PImat[,2]==pTime[j])
          if (length(aa)!=1){
            bb <- aa[abs(PImat[aa,3]-rMz)!=min(abs(PImat[aa,3]-rMz))]
            if(length(setdiff(aa,bb))!=1){
              cc <- setdiff(aa,bb)
              cc <- cc[PImat[cc,4]!=max(PImat[cc,4])]
              bb <- union(bb,cc)
            }
            PImat <- PImat[-bb,]
          }
        }
      }

      tt <- which(diff(PImat[,2])>5*tInl)
      if (length(tt>0)){
        uu <- union(1,tt+1)
        vv <- union(tt,nrow(PImat))
        ww <- which(vv-uu==max(vv-uu))[1]
        PImat <- PImat[uu[ww]:vv[ww],]
      }
      
      if(is.null(nrow(PImat))==T){ssid <- PImat[1]}else{ssid <- PImat[,1]}
      id <- setdiff(id[-1],ssid)
      
      if (length(ssid)<5|max(PImat[,4])<level){next}
      PIC <- unique(PImat[,c(2,4)])
      colnames(PIC) <- c("rt","intensity")
      rownames(PIC) <- NULL
      mzmax <- max(PImat[,3])
      mzmin <- min(PImat[,3])
      mz <- median(PImat[,3])
      rtmin <- min(PImat[,2])
      rtmax <- max(PImat[,2])
      maxinte <-max(PImat[,4])
      rt <- PImat[which(PImat[,4]==maxinte)[1],2]
      Time <- seq(from=min(PIC[,1]),to=max(PIC[,1]),by=tInl)
      
      PICs[[length(PICs)+1]] <- PIC
      Info <- rbind(Info,c(mz,mzmin,mzmax,rt,rtmin,rtmax,maxinte))
      colnames(Info) <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","Intensity")
    }else{
      id <- setdiff(id[-1],ssid)
      next}
  }
  result <- list(PICs=PICs,Info=Info)
  
  return (result)
}

getPIC = function(filename,range,level,alpha=0.3,Iterate=4000,n=3)
{
  library(Rcpp)
  library(stats)
  library(parallel)
  library(foreach)
  subPIC <- subPIC
  
  options(warn = -1)
  
  data <- LoadData(filename)
  mat <- data.frame(data$Mat)
  mat <- mat[order(mat[,2]),]
  tInl <- data$Info$freq
  mzlist <- mat[,2]
  mzdiff <- diff(mzlist)
  Dens <- density(mzdiff,kernel="gaussian",n=length(mzdiff)/5000,na.rm = T)
  reff <- mean(Dens$y)+sd(Dens$y)
  num <- which(Dens$y>reff)
  gap <- Dens$x[num[length(num)]]
  if (gap<=0){gap <- 0.005}
  
  library(doParallel)
  cl <- makeCluster(getOption("cl.cores", detectCores(logical = F)))
  registerDoParallel(cl)

  temp2 <- c(which(mzdiff>n*gap),length(mzdiff+1))
  temp1 <- c(1,temp2[1:length(temp2)-1]+1)
  nfrag <- length(temp1)
  
  datalist <- list()
  for (i in 1:nfrag){
    if (temp2[i]-temp1[i]>5){
      frag <- mat[temp1[i]:temp2[i],]
      if(max(frag[,3]>level)){datalist[[length(datalist)+1]] <- frag
      }
    }
  }
  
  PIC <-foreach(i=1:length(datalist)) %dopar%
    subPIC(datalist[[i]],range=range,level=level,alpha=alpha,Iterate=Iterate,tInl=tInl,gap=gap)
  stopCluster(cl)
  PIClist <- PIC[[1]]$PICs
  Infolist <- PIC[[1]]$Info
  for (kk in 2:length(PIC)){
    PIClist <- c(PIClist,PIC[[kk]]$PICs)
    Infolist <- rbind(Infolist,PIC[[kk]]$Info)
  }
  
  output <- list()
  output$PICs <- PIClist
  output$Info <- Infolist
  DataInfo <- data$Info
  output$rt <- seq(DataInfo$dStartTime,DataInfo$dEndTime,DataInfo$freq)
  return(output)
}

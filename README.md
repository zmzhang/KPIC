# KPIC  
  Pure Ion Chromatograms Extraction via Optimal K-means Clustering  

## Installation  

####Install Depends:  

    Rcpp, stats, parallel, Ckmeans.1d.dp, mzR, foreach, doParallel, iterators,
    devtools， Matrix

####Install KPIC:  

    library(devtools);  
    httr::set_config(httr::config(ssl_verifypeer = 0L))
    install_github("hcji/KPIC")

## Usage 

#### processing one sample
    path = 'E:/simulated_data.mzXML'
    PICs = getPIC(path,15,400,alpha=0.1)
    result = PICrefine(PICs,n=1)
    str(result)
    PICplot(result,1) 
	
#### processing a set of samples

    path = 'E:/samples'
	xset = KPICset(path,path,15,400,alpha=0.1,refine=1)
	xset = rtcor(xset,mztol=0.05,pw=20,lambda=1.5)
	str(xst)
  
## Contact
  ji.hongchao@foxmail.com
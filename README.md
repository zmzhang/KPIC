# KPIC  
  Pure Ion Chromatograms Extraction via Optimal K-means Clustering  

## Installation  

####Install Depends:  
	
    install.packages(c("devtools","Rcpp","Ckmeans.1d.dp","foreach","doParallel","iterators"))
	source("https://bioconductor.org/biocLite.R")
	biocLite("mzR")
	

####Install KPIC:  
	
    library(devtools)
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
	xset = KPICset(path,15,400,alpha=0.1,refine=1)
	xset = rtcor(xset,path)
	str(xset)
	
## Contact
  ji.hongchao@foxmail.com
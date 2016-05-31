PICplot = function(PICrefine,index)
{
  PIC <- PICrefine$PICs[[index]]
  rt <- PIC[,1]
  inte <- PIC[,2]
  plot(rt,inte,wd=2,pch=15,col='blue',cex.lab=1.5,cex.axis=1.3,font=2,
       type="o",xlab="Retention time",ylab="Intensity",
       main="Pure Ion Chromatography",cex.main=1.5)
}
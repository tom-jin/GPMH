####Summary for the class####
summary.GPMH=function(x){
  x=as.numeric(x)
  describe(x)[c(2:5,8:12)]
}
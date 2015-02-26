plot.GPMH <- function(object, ...){
  x <- obj$x
  xn <- data.frame(x=(x-mean(x))/sd(x))
  Xculm <- array(data=NA,dim=length(x))
  Xculm[1] <- x[1]
  
  for(i in 2:length(x)){
    Xculm[i] <- (x[i] + (i-1)*Xculm[i-1])/i
  }
  
  xculm <- data.frame(x=Xculm, simulation=1:length(x))
  p1 <- ggplot(xn, aes(x=x)) +
    geom_histogram(aes(y = ..density.., fill=..count..), binwidth = 0.15) +
    geom_density(fill=NA) +
    ggtitle("Histogram of MCMC Samples")
  
  p2<-ggplot(xculm, aes(simulation,x)) + geom_line() +
    ggtitle("MCMC Cumulative Means")
  
  grid.arrange(p1,p2,ncol=1)
}
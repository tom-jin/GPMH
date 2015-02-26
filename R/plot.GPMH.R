#' @import ggplot2
#' @export
plot.GPMH <- function(object, ...){
  
  #Source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  x <- object$x
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
  
  multiplot(p1, p2, cols=1)
}
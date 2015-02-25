######MH######
xMH<-MH(target=function(x) {0.2*dnorm(x,12,1)+0.3*dnorm(x,8,1)+0.4*dnorm(x,4,1)+0.6*dnorm(x,15,1)+0.3*dnorm(x,20,1)},
         kernel=function(x) {rnorm(n=1, x)},
         dkernel=function(x,y) {dnorm(y, x)},
         init.state=5,
         n=100000)

plot(xMH$x)
summary(xMH$x)

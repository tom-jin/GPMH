#Generalised MH algorithm WITH SOME C CALLING from Cgmh.c

rcGPMH<-function (target,
                 kernel,  # to draw from the distribution of mean = x
                 dkernel, # to evaluate in y the distribution of mean = x 
                 init.state,
                 n,
                 N=8)  # points to "parallely" draw from the proposal at each iteration
  
{
  d <- length(init.state)
  I <- 1 
  X <- matrix(NA, nrow = n*(N + 1), ncol = d)
  X[1, ] <- init.state
  Y <- matrix(NA, nrow = N+1, ncol = d)
  Y[1,] <- X[1,]
  for(j in 2:(N+1)){
    Y[j,]<-kernel(X[1, ]) # generate N new points from the proposal
  }
  
  A<-numeric(N+1) 
  
  for (i in 1:n) { 
    
    A<-cgmh(A,Y,N)  # Call the C chunk to compute the vector of probabilities A
    
    # sample from the N+1 points to obtain the new N+1 MCMC samples 
    
    X[((i-1)*N+i):(i*(N + 1)),]<-Y[sample(seq(1:(N+1)),replace=TRUE,prob=A),]
    
    # Now we set Y for the next iteration:
    I<-sample(((i-1)*N+i):(i*(N + 1)),1)  # draw uniformly the index I...
    Y[1,]<-X[I, ] # ...of the X[I,] that will generate N new points
    for(j in 2:(N+1)){
      Y[j,]<-kernel(X[I, ]) # generate N new points
    }
  }
  
  output <- list(x = X[-1, ], n = n, call = match.call())
  class(output) <- "GPMH"
  return(output)
}
#  compare 3 approaches to solving a linear system of ODEs
#  METHOD 1) Numerical solver
#  METHOD 2) Matrix exponential formula implemented in R 
#  METHOD 3) Derived analytic solution
# All approaches give the same result.

rm(list=ls())


library(deSolve)
library(tidyr)
library(ggplot2)
library(expm)

## set up the A matrix
W <- matrix(c(-.2,0,.3,-.5),nrow=2,byrow=TRUE)
W

diag(nrow(W))-W
solve(diag(nrow(W))-W)
eigen(diag(nrow(W))-W)

B <- matrix(c(2,0,1,
              0,1,1),nrow=2,byrow=TRUE)
temp <- matrix(0,nrow=3,ncol=5)
A <- rbind(temp,cbind(B,W))


### METHOD 1)

## set up the necessary functions for deSolve
LinTrans <- function(t, state, parameters) {
  return(list(colSums(t(parameters)*state)))
}

state <- c(1,rep(0,nrow(A)-1))
times <- seq(0, 100, by = 0.01)


LinTrans(7,state,A)
LinTrans(10,state,A)
A%*%matrix(state,ncol=1)


out <- ode(y=state,times=times,func=LinTrans,parms =A)
head(out)


## plot deSolve output
df <- pivot_longer(as.data.frame(out),2:ncol(out))
ggplot(df,aes(x=time,y=value,col=name)) +
  geom_line()
  
### METHOD 2)
## compute limiting distribution of system
t <- 100
expm(t*A)

### METHOD 3)
## analytic solution for lower left block element
solve(-W)%*%B

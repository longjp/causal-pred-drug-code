### compare a) causal discovery vs b) regression approach
### compare in 3 cases
### SIM 1) Random train-test split (both a,b are good)
### SIM 2) Random train-test split with misspecified B matrix (b is good)
### SIM 3) Leave one drug out (a is good)
rm(list=ls())
library(ggplot2)
library(igraph)
set.seed(1234)

### plotting parameteres
## margin size for plots
mar_size <- 0
## plot width
ll <- 1.2
## threshold for plotting edge
thres <- 0.2

rr <- 5
a <- t(combn(rr,2))
b <- matrix(0,nrow=nrow(a),ncol=rr)
a <- rbind(cbind(1:nrow(a),a[,1]),cbind(1:nrow(a),a[,2]))
b[a] <- 1
Bt <- rbind(diag(rr),b/2)
B <- t(Bt)


## use all combinations of drugs
a <- t(combn(ncol(B),2))
b <- matrix(0,nrow=nrow(a),ncol=ncol(B))
a <- rbind(cbind(1:nrow(a),a[,1]),cbind(1:nrow(a),a[,2]))
b[a] <- 1
D <- b
D <- D[sample(nrow(D)),]


A <- matrix(c(0,0,0,0,0,
              1.6,0,0,0,0,
              1.2,0,0,0,0,
              0,0,2,0,0,
              0,0,0,0,0),nrow=5,byrow=TRUE)

W <- A - diag(nrow(A))


# build the graph object
colnames(A) = rownames(A) = paste0("X",1:nrow(A))
network <- graph_from_adjacency_matrix(t(A),weighted=TRUE)
coords <- layout_in_circle(network)

# plot true network

pdf("figs/network_true.pdf")
par(mar=rep(mar_size,4))
plot(network,layout=coords,
     edge.width=3*E(network)$weight,
     edge.arrow.size=1.5,
     edge.arrow.width=1.5,
     vertex.size=40,
     vertex.label.cex=2,
     xlim=c(-ll,ll),
     ylim=c(-ll,ll),
     edge.curved=TRUE)
dev.off()








# W <- matrix(c(-.2,0,0,0,0,
#               .3,-.5,0,0,0,
#               0.2,0,-.7,0,0,
#               0,0,0.5,-0.2,0,
#               0,0,0,0,-.1),nrow=5,byrow=TRUE)
solve(diag(nrow(W))-W)



X <- D%*%t(B)%*%(-t(solve(W)))
X <- X + matrix(rnorm(prod(dim(X)),sd=.2),nrow=nrow(X))


temp <- D%*%t(B)
Tt <- solve(t(temp)%*%temp)%*%t(temp)%*%X
What <- -solve(t(Tt))
sum(abs(What-W))





## SIM 1) random training / test set split

ntr <- floor(2*nrow(X)/3)
## training / test set splits
Xtr <- X[1:ntr,]
Dtr <- D[1:ntr,]
Xte <- X[(ntr+1):nrow(X),]
Dte <- D[(ntr+1):nrow(D),]

## predict using regression
R <- solve(t(Dtr)%*%Dtr)%*%t(Dtr)%*%Xtr
Xpr <- Dte%*%R

## predict using cellbox
temp <- Dtr%*%t(B)
Tt <- solve(t(temp)%*%temp)%*%t(temp)%*%Xtr
What <- -solve(t(Tt))
Ahat <- diag(nrow(What)) + What
sum(abs(What-W))
Xpr2 <- Dte%*%t(B)%*%Tt




# plot estimate A
Ahatp <- Ahat
Ahatp[abs(Ahatp)<thres] <- 0

colnames(Ahatp) = rownames(Ahatp) = paste0("X",1:nrow(Ahatp))
network <- graph_from_adjacency_matrix(t(Ahatp),weighted=TRUE)
coords <- layout_in_circle(network)

# plot true network

pdf("figs/network_RF.pdf")
par(mar=rep(mar_size,4))
plot(network,layout=coords,
     edge.width=abs(3*E(network)$weight),
     edge.arrow.size=1.5,
     edge.arrow.width=1.5,
     vertex.size=40,
     vertex.label.cex=2,
     xlim=c(-ll,ll),
     ylim=c(-ll,ll),
     edge.curved=TRUE)
dev.off()






## quality of result
pdf("figs/sim1-1.pdf",width=8,height=4)
par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
plot(as.vector(Xte),as.vector(Xpr),xlab="True",ylab="Regression Prediction")
plot(as.vector(Xte),as.vector(Xpr2),xlab="True",ylab="Causal Prediction")
dev.off()


cor(as.vector(Xte),as.vector(Xpr2))
cor(as.vector(Xte),as.vector(Xpr))

A
Ahat*(abs(Ahat)>0.1)

## assemble into single data frame
out <- data.frame(est=c(rep("Regression",prod(dim(Xte))),rep("Causal",prod(dim(Xte)))),
                  exp="RF",
                  true=c(as.vector(Xte),as.vector(Xte)),
                  pred=c(as.vector(Xpr),as.vector(Xpr2)))




## SIM 2) 
## misspecified B
## training / test set splits
ntr <- floor(2*nrow(X)/3)
## training / test set splits
Xtr <- X[1:ntr,]
Dtr <- D[1:ntr,]
Xte <- X[(ntr+1):nrow(X),]
Dte <- D[(ntr+1):nrow(D),]

## predict using regression
R <- solve(t(Dtr)%*%Dtr)%*%t(Dtr)%*%Xtr
Xpr <- Dte%*%R

## predict using cellbox

## misspecified B
Bm <- t(do.call(rbind,replicate(3,diag(rr),simplify=FALSE)))
temp <- Dtr%*%t(Bm)
Tt <- solve(t(temp)%*%temp)%*%t(temp)%*%Xtr
What <- -solve(t(Tt))
Ahat <- diag(nrow(What)) + What
sum(abs(What-W))
Xpr2 <- Dte%*%t(Bm)%*%Tt



# plot estimate A
Ahatp <- Ahat
Ahatp[abs(Ahatp)<thres] <- 0

colnames(Ahatp) = rownames(Ahatp) = paste0("X",1:nrow(Ahatp))
network <- graph_from_adjacency_matrix(t(Ahatp),weighted=TRUE)
coords <- layout_in_circle(network)

# plot true network

pdf("figs/network_RF_missB.pdf")
par(mar=rep(mar_size,4))
plot(network,layout=coords,
     edge.width=abs(3*E(network)$weight),
     edge.arrow.size=1.5,
     edge.arrow.width=1.5,
     vertex.size=40,
     vertex.label.cex=2,
     xlim=c(-ll,ll),
     ylim=c(-ll,ll),
     edge.curved=TRUE)
dev.off()







## quality of result
pdf("figs/sim1-2.pdf",width=8,height=4)
par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
plot(as.vector(Xte),as.vector(Xpr),xlab="True",ylab="Regression Prediction")
plot(as.vector(Xte),as.vector(Xpr2),xlab="True",ylab="Causal Prediction")
dev.off()

cor(as.vector(Xte),as.vector(Xpr2))
cor(as.vector(Xte),as.vector(Xpr))



## assemble into single data frame
out_temp <- data.frame(est=c(rep("Regression",prod(dim(Xte))),rep("Causal",prod(dim(Xte)))),
                      exp="RF with B Misspecified",
                      true=c(as.vector(Xte),as.vector(Xte)),
                      pred=c(as.vector(Xpr),as.vector(Xpr2)))

out <- rbind(out,out_temp)




## SIM 3)
## leave one drug out CV
Xpr_cell <- matrix(NA_real_,ncol=ncol(X),nrow=nrow(X))
Xpr_reg <- Xpr_cell
for(ii in 1:ncol(D)){
  ix <- D[,ii]!=0
  
  Xtr <- X[!ix,]
  Dtr <- D[!ix,]
  Dte <- D[ix,]
  
  ## predict using regression
  R <- lm.fit(Dtr,Xtr)$coefficients
  R[is.na(R)] <- 0
  Xpr_reg[ix,] <- Dte%*%R
  
  ## predict using cellbox
  temp <- Dtr%*%t(B)
  Tt <- solve(t(temp)%*%temp)%*%t(temp)%*%Xtr
  What <- -solve(t(Tt))
  Ahat <- diag(nrow(What)) + What
  sum(abs(What-W))
  Xpr_cell[ix,] <- Dte%*%t(B)%*%Tt
}
  



# plot estimate A
Ahatp <- Ahat
Ahatp[abs(Ahatp)<thres] <- 0

colnames(Ahatp) = rownames(Ahatp) = paste0("X",1:nrow(Ahatp))
network <- graph_from_adjacency_matrix(t(Ahatp),weighted=TRUE)
coords <- layout_in_circle(network)

# plot true network

pdf("figs/network_LODO.pdf")
par(mar=rep(mar_size,4))
plot(network,layout=coords,
     edge.width=abs(3*E(network)$weight),
     edge.arrow.size=1.5,
     edge.arrow.width=1.5,
     vertex.size=40,
     vertex.label.cex=2,
     xlim=c(-ll,ll),
     ylim=c(-ll,ll),
     edge.curved=TRUE)
dev.off()





## quality of result
pdf("figs/sim1-3.pdf",width=8,height=4)
par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
plot(as.vector(X),as.vector(Xpr_reg),xlab="True",ylab="Regression Prediction")
plot(as.vector(X),as.vector(Xpr_cell),xlab="True",ylab="Causal Prediction")
dev.off()

## assemble into single data frame
out_temp <- data.frame(est=c(rep("Regression",prod(dim(X))),rep("Causal",prod(dim(X)))),
                      exp="LODO",
                      true=c(as.vector(X),as.vector(X)),
                      pred=c(as.vector(Xpr_reg),as.vector(Xpr_cell)))

out <- rbind(out,out_temp)


## order factor levels
out$est <- factor(out$est,levels=c("Regression","Causal"))
out$exp <- factor(out$exp,levels=c("RF","RF with B Misspecified","LODO"))

## make nice plot
pp <- ggplot(out, aes(x = true, y = pred)) + geom_point(alpha=0.4)
pp <- pp + facet_grid(est ~ exp) + xlab("True Response") + ylab("Predicted Response")

pdf("figs/sim1-4.pdf",width=8,height=4.3)
print(pp)
dev.off()


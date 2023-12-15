### compare a) regression and b) CSR/Cellbox
### SIM 1) Random train-test split (both a,b are good)
### SIM 2) Random train-test split with misspecified B matrix (a is good)
### SIM 3) Leave one drug out (b is good)
rm(list=ls())
library(ggplot2)
library(dplyr)
library(igraph)
set.seed(20230607)


outf <- "sim-dag-figs/"
unlink(outf,recursive=TRUE)
dir.create(outf,recursive=TRUE)


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


# plot true network
# build the graph object
colnames(A) = rownames(A) = paste0("X",1:nrow(A))
network <- graph_from_adjacency_matrix(t(A),weighted=TRUE)
coords <- layout_in_circle(network)
pdf(paste0(outf,"network_true.pdf"))
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


## generate data
ep <- matrix(rnorm(nrow(D)*nrow(B),sd=.1),nrow=nrow(D))
X <- (D%*%t(B) + ep)%*%t(solve(diag(nrow(A))-A))


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
Ahat <- t(diag(nrow(Tt))-solve(Tt))
Xpr2 <- Dte%*%t(B)%*%Tt

# plot estimate A
Ahatp <- Ahat
Ahatp[abs(Ahatp)<thres] <- 0
colnames(Ahatp) = rownames(Ahatp) = paste0("X",1:nrow(Ahatp))
network <- graph_from_adjacency_matrix(t(Ahatp),weighted=TRUE)
coords <- layout_in_circle(network)
pdf(paste0(outf,"network_RF.pdf"))
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
pdf(paste0(outf,"RF.pdf"),width=8,height=4)
par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
plot(as.vector(Xte),as.vector(Xpr),xlab="True",ylab="Regression Prediction")
plot(as.vector(Xte),as.vector(Xpr2),xlab="True",ylab="Causal Prediction")
dev.off()


#cor(as.vector(Xte),as.vector(Xpr2))
#cor(as.vector(Xte),as.vector(Xpr))


## assemble into single data frame
out <- data.frame(est=c(rep("LR",prod(dim(Xte))),rep("CSR",prod(dim(Xte)))),
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
Ahat <- t(diag(nrow(Tt))-solve(Tt))
Xpr2 <- Dte%*%t(Bm)%*%Tt

# plot estimate A
Ahatp <- Ahat
Ahatp[abs(Ahatp)<thres] <- 0

colnames(Ahatp) = rownames(Ahatp) = paste0("X",1:nrow(Ahatp))
network <- graph_from_adjacency_matrix(t(Ahatp),weighted=TRUE)
coords <- layout_in_circle(network)
pdf(paste0(outf,"network_RF_missB.pdf"))
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
pdf(paste0(outf,"RF-Bmiss.pdf"),width=8,height=4)
par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
plot(as.vector(Xte),as.vector(Xpr),xlab="True",ylab="Regression Prediction")
plot(as.vector(Xte),as.vector(Xpr2),xlab="True",ylab="Causal Prediction")
dev.off()

cor(as.vector(Xte),as.vector(Xpr2))
cor(as.vector(Xte),as.vector(Xpr))



## assemble into single data frame
out_temp <- data.frame(est=c(rep("LR",prod(dim(Xte))),rep("CSR",prod(dim(Xte)))),
                      exp="RF with B Misspecified",
                      true=c(as.vector(Xte),as.vector(Xte)),
                      pred=c(as.vector(Xpr),as.vector(Xpr2)))

out <- rbind(out,out_temp)




## SIM 3)
## leave one drug out CV
ntr <- floor(2*nrow(X)/3)
## training / test set splits
Xtr <- X[1:ntr,]
Dtr <- D[1:ntr,]
Xte <- X[(ntr+1):nrow(X),]
Dte <- D[(ntr+1):nrow(D),]


Xpr_reg <- matrix(NA_real_,nrow=nrow(Xte),ncol=ncol(Xte))
Xpr_cell <- matrix(NA_real_,nrow=nrow(Xte),ncol=ncol(Xte))
for(ii in 1:nrow(Xte)){
  ## select 1 test observation
  Dte1 <- Dte[ii,,drop=FALSE]
  ## choose one of the two drug used on test obs (at random)
  LOD <- sample(which(Dte1==1),1)
  ## find training where left out drug is 0
  ix <- Dtr[,LOD]==0
  ## created training set where LOD not used
  Xtr_sub <- Xtr[ix,]
  Dtr_sub <- Dtr[ix,]
  ## predict using regression
  R <- lm.fit(Dtr_sub,Xtr_sub)$coefficients
  R[is.na(R)] <- 0
  Xpr_reg[ii,] <- Dte1%*%R
  ## predict using cellbox
  temp <- Dtr_sub%*%t(B)
  Tt <- solve(t(temp)%*%temp)%*%t(temp)%*%Xtr_sub
  Ahat <- t(diag(nrow(Tt))-solve(Tt))
  Xpr_cell[ii,] <- Dte1%*%t(B)%*%Tt
}
  



# plot estimate A
Ahatp <- Ahat
Ahatp[abs(Ahatp)<thres] <- 0
colnames(Ahatp) = rownames(Ahatp) = paste0("X",1:nrow(Ahatp))
network <- graph_from_adjacency_matrix(t(Ahatp),weighted=TRUE)
coords <- layout_in_circle(network)
pdf(paste0(outf,"network_LODO.pdf"))
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
pdf(paste0(outf,"LODO.pdf"),width=8,height=4)
par(mfcol=c(1,2),mar=c(4.5,4.5,1,1))
plot(as.vector(Xte),as.vector(Xpr_reg),xlab="True",ylab="Regression Prediction")
plot(as.vector(Xte),as.vector(Xpr_cell),xlab="True",ylab="Causal Prediction")
dev.off()

## assemble into single data frame
out_temp <- data.frame(est=c(rep("LR",prod(dim(Xte))),rep("CSR",prod(dim(Xte)))),
                      exp="LODO",
                      true=c(as.vector(Xte),as.vector(Xte)),
                      pred=c(as.vector(Xpr_reg),as.vector(Xpr_cell)))

out <- rbind(out,out_temp)


## order factor levels
out$est <- factor(out$est,levels=c("LR","CSR"))
out$exp <- factor(out$exp,levels=c("RF","RF with B Misspecified","LODO"))

## make nice plot
pp <- ggplot(out, aes(x = true, y = pred)) + geom_point(alpha=0.4)
pp <- pp + facet_grid(est ~ exp) + xlab("True Response") + ylab("Predicted Response") +
  geom_abline(slope=1,intercept=0,alpha=0.5)

rho2 <- out %>%
  group_by(est,exp) %>%
  summarise(rho = round(cor(true,pred), 3))


mae2 <- out %>%
  group_by(est,exp) %>%
  summarise(mae = round(mean(abs(true-pred)),3))


pp <- pp + geom_text(x = 1, y = 3.5, 
            aes(label = paste0("MAE: ", mae)), 
            data = mae2)
pp <- pp + geom_text(x = 1, y = 4, 
            aes(label = paste0("Pearson: ", rho)), 
            data = rho2)

pdf(paste0(outf,"all.pdf"),width=8,height=4.3)
print(pp)
dev.off()



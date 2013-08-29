# Colbert Seanker HW 5 PCA and Factor Analysis
A <-  as.matrix(read.table("population_a.txt", sep=" "))
B <-  as.matrix(read.table("population_b.txt", sep=" "))
# returns normalized covariance matrix
normalize <- function(C) {
  mu <- apply(C, 2, sum)/nrow(C); p <- mu/2
  M  <- C-mu/(p*(1-p))^(.5)
  return(list(norm=M %*% t(M), # normalized and mean centered
              center=(C-mu) %*% t(C-mu) # mean centered, not normalized
             )
        ) 
}
A.norm <- normalize(A)$norm; B.norm <- normalize(B)$norm; 
pca.plot <- function(x, group) {
  pc  <- list(one=eigen(x)$vectors[,1], two=eigen(x)$vectors[,2]) # list of principal components
  lim <- c(pc$one, pc$two)
  plot(pc$one, type="p", col="red",xlab="",ylab="",ylim=range(lim)); par(new=TRUE);
  plot(pc$two, type="p", col="blue",xlab="individual",ylab="loading value",
  main=paste("Population", group, "PCA"));  
  return(pc)
} 
png(file="PCA plots.png", width=512, height=512,res=72)
op <- par(mfrow=c(2,1));
A.pc <- pca.plot(A.norm,"A");  B.pc <-pca.plot(B.norm,"B");
dev.off()

# Compute factors
A.c  <- normalize(A)$center; B.c <- normalize(B)$center
A.f2 <- factanal(cov=A.c,factors=2)$loadings 
B.f2 <- factanal(cov=B.c,factors=2)$loadings
plot.fac  <- function(factors, group) {#plot factor A vs factor B
     plot(factors[,1], factors[,2], type="p",xlab="factor 1",ylab="factor 2",
     main=paste("factor 1 vs factor 2", group));
}
png(file="Factor1 vs Factor2 plots.png", width=512, height=512,res=72)
op <- par(mfrow=c(2,1));
plot.fac(A.f2,"A"); plot.fac(B.f2,"B")
dev.off()
# Plot factors loadings
png(file="Factor Analysis plots.png", width=512, height=512,res=72)
op <- par(mfrow=c(3,2));
A.f3 <- list(fac=factanal(cov=A.c,factors=3)$loadings, group="A")
B.f3 <- list(fac=factanal(cov=B.c,factors=3)$loadings, group="B")
for (facs in list(A.f3, B.f3)){
    for(i in 1:3){     
      plot(facs$fac[,i], type="p", xlab="individual", ylab="loading value",
           main=paste("Population", facs$group, "factor", i));
    }
}
dev.off()

# The factors reveal much more about population structure

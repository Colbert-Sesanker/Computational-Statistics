# EM contour plots by Colbert Sesanker March 2013, source this file for a demo
# EM demo on 2d normals, assuming 2 latent variables, plots the gaussian contours that fit data 

library(mvtnorm)
pts_mat <- as.matrix(read.table("points_hw4.txt", sep=" "))
points  <- split(pts_mat, c(row(pts_mat))) # split matrix into a mapable list of points
mc_pts  <- function(mean) Map(function(x) x-mean, points)   # mean centered points split into mappable list
n     <- length(points)
nth   <- function(n, tuple_list) Map(function(x) x[n], tuple_list) # returns list of the nth element of each tuple from list of tuples
mat.sum <- function(mat, weights=rep(1,length(mat))) Reduce("+", Map("*", mat, weights)) # weighted matrix sum of a list of matricies/vectors
theta <- list(
             pi= c(0.5,0.5),
             mu1=c(2.5,7),
             mu2=c(3.6,12),
             sigma1=matrix(c(3,5,5,15),ncol=2),
             sigma2=matrix(c(3,5,5,15),ncol=2)
             )

EM <- function(theta) {
  E.step <- Map(function(x,y) c(x,y)/(x+y),
              theta$pi[1] * dmvnorm(pts_mat,mean=theta$mu1,sigma=theta$sigma1),
              theta$pi[2] * dmvnorm(pts_mat,mean=theta$mu2,sigma=theta$sigma2)
               )
  z.eq1 <-nth(1, E.step); z.eq0 <-nth(2, E.step) # Posterior probabilities of z = 0 or 1
  z1 <- mat.sum(z.eq1);   z0 <- mat.sum(z.eq0); 
  mu1=mat.sum(points, z.eq1)/z1; mu2=mat.sum(points, z.eq0)/z0
  M.step <- list(
                pi= c(z1/n, z0/n),
                mu1=mu1,
                mu2=mu2,
                sigma1=mat.sum(Map(function(x) x %o% x, mc_pts(mu1)), z.eq1)/z1, # This follows from the definition of sigma using the outer product
                sigma2=mat.sum(Map(function(x) x %o% x, mc_pts(mu2)), z.eq0)/z0 
                )
  return(M.step)
}

contour.plot <- function(theta, iter) {
  x=y=seq(1,6.5,.01)
  mixture.pdf <- function(x,y) theta$pi[1] * dmvnorm(cbind(x,y),mean=theta$mu1,sigma=theta$sigma1) +
                               theta$pi[2] * dmvnorm(cbind(x,y),mean=theta$mu2,sigma=theta$sigma2) 
  mixture.contour <- outer(x,y, mixture.pdf)
  contour(x, y, mixture.contour, nlevels=10, drawlabel=FALSE, col="blue",
          main=paste("EM estimation",iter, "iterations"),xlab="x", ylab="y")
          points(pts_mat)
}
for (iter in 1:50){
 theta <- EM(theta) 
 contour.plot(theta, iter) 
}
png(file="EM_mixture_contour_plot", width=512, height=512,res=72)
contour.plot(theta, 50) 
dev.off()
















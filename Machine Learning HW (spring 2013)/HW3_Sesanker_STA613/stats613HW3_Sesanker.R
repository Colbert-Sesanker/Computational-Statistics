### STA 613 HW 3, Colbert Sesanker

library(MASS) # for ginv

genes <- as.matrix(read.table("hw1_genes.txt", sep=" "))
SNPs <- as.matrix(read.table("hw1_genotypes.txt", sep=" "))

# Question 2a

# This returns the mean or the posterior mean, beta_tilde = inv((XX' + (1/lambda)*eye(#covariates))((mu/lambda)*eye(#covariates) + X'y)
# X is the design matrix, y the response, l is the prior variance on beta, mu is the prior mean, likelihood variance assumed = 1
bay_regress <- function(X, y, l, mu) return( ginv(t(X)%*%X + (1/l)*diag(5)) %*% ((1/l)*diag(5)%*%rep(mu,5) + t(X)%*%y)) # morse penrose inverse with prior weight

freq_regress <- function(X, y) ginv(X)%*%y # the beta_hat is the morse penrose inverse of design matrix X (MLE estimate) times reponse, y
                                           # only difference from above baysian estimate is the diag matrix in the pseudo inverse. note: the ginv in  
                                           # bay_regress is just the standard non-generalized inverse 
bay_models <- array(0, dim = c(5,10))
freq_models <- array(0, dim = c(5,10))
for (i in 1:10) {
  bay_models[,i] <- bay_regress(SNPs, genes[,i], 1, 100) # tune mu and lambda parameters here
  freq_models[,i] <- freq_regress(SNPs, genes[,i])
}

# Plots 2b
png(file=" Baysian gene-SNPs regression mu = 100",width=1000,height=1000,res=72)
op <- par(mfrow=c(5,2)); 
par(mar = rep(2, 4))
for (i in 1:10) {
  plot(SNPs %*% bay_models[ ,i], genes[,i],  main=paste("gene", i)); abline(0,1)
}
dev.off()

png(file=" Frequentist gene-SNPs regression",width=1000,height=1000,res=72)
op <- par(mfrow=c(5,2)); 
for (i in 1:10) {
  # plot(SNPs %*% freq_models[ ,i], genes[,i])
  plot(SNPs %*% freq_models[ ,i], genes[,i],  main=paste("gene", i)); abline(0,1)
}
dev.off()
# Question 2b
# The baysian estimates are significantly more dispersed along the vertical axis genes. The bands are also narrowersalong the SNPs axis in the baysian # # # regression.
# 2c: increasing lambda takes away the impact of the prior, this is cleat from the equations as the posterior mean becomes the MLE as lambda -> inf
#     decreasing lambda increases the weight of the prior and narrows the estimation bands
# 2d: increasing mu narrows the bands along the SNPs axis, similar, but less significant, effect to decreasing lambda

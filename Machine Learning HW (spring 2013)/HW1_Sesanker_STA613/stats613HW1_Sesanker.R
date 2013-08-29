## HW 1
## Colbert Sesanker, CBB 540

### Question 1
## 1 a
n <- 1000
samples_p <- rpois(n, 23)
hist(samples_p)
above <- function(X) 
{
if (50 < X) return(TRUE) else return(FALSE)
}
below <- function(X) 
{
if (X < 10) return(TRUE) else return(FALSE)
}
length(Filter(above,samples_p)) # got 0
length(Filter(below,samples_p)) # got 2
var(samples_p) # got 21.24

## 1 b
n <- 1000
samples_n <- rnbinom(n, 23, .5)
hist(samples_n)
length(Filter(above, samples_n)) # got 1
length(Filter(below, samples_n)) # got 8, bigger than pois
var(samples_n) # got 45.23, bigger than pois
## 1 c
sort_pois <- sort(samples_p)
sort_neg <- sort(samples_n)
plot(sort_pois)
"
poission makes a regular sideways S curve
"
plot(sort_neg)
"
same shape as pois, but takes on large and small
values more often

"
## Question 2 on paper

### Question 3
## a
genes <- as.matrix(read.table("hw1_genes.txt", sep=" "))
num_genes <- dim(genes)[2]
op <- par(mfrow=c(num_genes,1))
par(mar = rep(2, 4)) # resize margins to fit plots
for (i in 1:num_genes) {
  gene <- genes[ ,i]
  hist(gene, main=paste("gene", i, "expression in 500 individuals" ))
}
"
genes 1,2,6,7 look off from normal, 4,5 look close to normal
"
## 3 b
# Calculate means and variances for each gene:
mean <- c(); sd <- c(); 
for (i in 1:num_genes) {
  gene <- genes[ ,i]
  mean[i] <- mean(gene)
  sd[i] <- sd(gene)
}
# draw 10 samples of 500 from normal
norm_sample <- array(0, dim = c(10,500))
for (i in 1:10) {
  norm_sample[i, ] <- rnorm(500, mean[i], sd[i])
}  
# qqplot
op <- par(mfrow=c(num_genes,1))
par(mar = rep(2, 4))
for (i in 1:num_genes) {
  qqplot(norm_sample[i, ], genes[ ,i], main=paste("gene", i))
  abline(0,1)
}
"
qqplot yeilds same results as the intuitions from the histogram

"

## 3 c, standardize and compare 1 and 4
op <- par(mfrow=c(4,1))
for (i in c(1,4)) {
  st_gene <- (genes[ ,i] - mean[i])/ sd[i]
  qqplot(norm_sample[i, ], genes[ ,i], main=paste("gene", i)); abline(0,1)
  qqplot(rnorm(500, 0, 1), st_gene, main=paste("standardized gene", i)); abline(0,1)
}
"
standardization helped in some parts of gene one and made others worse
over all in seemed to allgin to the curve to y=x but did not significantly help outliers
"

## 3 d, standardize with qqnorm and compare 1 and 4
op <- par(mfrow=c(4,1))
for (i in c(1,4)) {
  std_gene <- qqnorm(genes[ ,i], plot.it = FALSE)
  std_gene <- unlist(std_gene[1])
  qqplot(norm_sample[i, ], genes[ ,i], main=paste("gene", i)); abline(0,1)
  qqplot(rnorm(500, mean(std_gene), sd(std_gene)), std_gene, main=paste("quantile projection  gene", i)); abline(0,1)
}

"
with quantile projection, the majority of the points lie on the line y=x save some outliers 
it does not help with extreme outliers
"
## 3 e covariance
cv <- cov(genes) #
myImagePlot(cv)
"
genes 6,6 | 5,5 | 1,3| and 10,10 covary / have high variance
"
## 3 f correlation
cr <- cor(genes) #
myImagePlot(cr)
"
genes  1,4 | 1,2 | 3,4 | 6,7 and 6,8 are correlated the most
"
### Question 4 a
SNPs <- as.matrix(read.table("hw1_genotypes.txt", sep=" "))
assoc <- array(0, dim = c(10,5))
for (i in 1:10) {
 for(j in 1:5) {
   s <- summary(lm(genes[ ,i] ~ SNPs[ ,j]))$fstatistic # extract p-value from model
   p_value <- pf(s[1],s[2],s[3],lower.tail=F)
   if (p_value <= .05) { # Check for significance at the 5% level
     assoc[i,j] = p_value # assoc matrix has the p-value for associated gene-snps and 0 otherwise
   }
 }
} 
# assoc matrix genes vs. SNPs association: SNPs on columns genes on rows
" 
 assoc
             [,1]         [,2]         [,3]        [,4]        [,5] 
 [1,] 0.000000000 1.858672e-09 1.901004e-03 0.000000000 0.000000000
 [2,] 0.000000000 4.372846e-02 0.000000e+00 0.000000000 0.000000000
 [3,] 0.000000000 2.031232e-22 2.317863e-03 0.000000000 0.000000000
 [4,] 0.000000000 8.432365e-07 4.073885e-10 0.000000000 0.000000000
 [5,] 0.007742018 0.000000e+00 0.000000e+00 0.000654598 0.000000000
 [6,] 0.000000000 0.000000e+00 0.000000e+00 0.000000000 0.000000000
 [7,] 0.000000000 0.000000e+00 0.000000e+00 0.000000000 0.000000000
 [8,] 0.000000000 0.000000e+00 0.000000e+00 0.000000000 0.000000000
 [9,] 0.000000000 0.000000e+00 0.000000e+00 0.000000000 0.000000000
[10,] 0.000000000 0.000000e+00 0.000000e+00 0.000000000 0.005109664

}

we see gene-SNPs (respectively): 1,2| 1,3| 2,2| 2,3| 3,2 | 3,3 | 4,2 |4,3| 5,1| 5,4 | 10,5 
are associated
"
## Question 4b
png(file="gene-SNP association",width=256,height=512,res=72)
op <- par(mfrow=c(2,1))
plot(SNPs[ ,3], genes[ ,4]); abline(lm(genes[ ,4] ~ SNPs[ ,3])) # associated
plot(SNPs[ , 1], genes[ ,1]); abline(lm(genes[ ,1] ~ SNPs[ ,1])) # not associated
dev.off()
## Question 4c
RSS1  <- summary(lm(genes[ ,4] ~ SNPs[ ,3]))$sigma # RSS: .667
RSS2  <- summary(lm(genes[ ,1] ~ SNPs[ ,1]))$sigma # RSS: 1.30195
r_sq1 <- summary(lm(genes[ ,4] ~ SNPs[ ,3]))$r.squared # r^2: .0755
r_sq2 <- summary(lm(genes[ ,1] ~ SNPs[ ,1]))$r.squared #r^2: .0015

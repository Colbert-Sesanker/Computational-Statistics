### STA 613 HW 2, Colbert Sesanker

genes <- as.matrix(read.table("hw1_genes.txt", sep=" "))
SNPs <- as.matrix(read.table("hw1_genotypes.txt", sep=" "))
phenotypes <- as.matrix(read.table("hw2_casecontrol.txt", sep=" "))

## Question 4 a, Computes sorted p_values for given genotype data
compute_p_vals <- function(SNPs) {
  p_vals <- array(0, dim = c(50))
  for (i in 1:10) {
    for (j in 1:5)  {
      s <- summary(lm(genes[ ,i] ~ SNPs[ ,j]))$fstatistic # extract p-value from model
      p_val <- pf(s[1],s[2],s[3],lower.tail=F)  
      p_vals[(j-1)*10 + i] <- p_val # array of p-values from 50 comparisons
     }
  }
  return(sort(p_vals)) # return p_values in sorted order
}

# p-values for non permuted genotypes
p_vals <- compute_p_vals(SNPs)

# returns True if p-value is less than threshold t in (0,1]
less_p_val <- function(t, p_val) {if (p_val <= t) return(TRUE) else return(FALSE)}
# returns number of p-values less than t
npv <- function(t, p_vals) length(Filter(partial(less_p_val, t), p_vals))

## Question 4c
npv(.001, p_vals) # returns 5, namely: 1.858672e-09, 2.031232e-22, 8.432365e-07, 4.073885e-10, 6.545980e-04

## Question 4d, permute genotype data by rows
permuted_SNPs <- SNPs[sample.int(nrow(SNPs)),]
permuted_p_vals <- compute_p_vals(permuted_SNPs)
# qq-plot against the actual p-values
png(file="actual vs null p-value association",width=512,height=512,res=72)
qqplot(p_vals, permuted_p_vals, main=paste("Actual vs null p-values")); abline(0,1)
dev.off()

## 4e, find the largest threshhold t s.t. #{permuted-p-vals <= t} / #{p-vals <= t} <= .1
find_threshold <- function(p_vals) {
  for (t in seq(1, .001, by= -.001)) {
  fdr <- npv(t, permuted_p_vals) / npv(t, p_vals) 
    if (fdr <= .1) {
       return(t)
    }
  }
 return(FALSE) # return false if no t satisfies
}
t <- find_threshold(p_vals) # 4f returns a threshhold of .011 for both permutations, (.06 on a third permutation)

# Question 5
compute_odds <- function() {
 odds_ratios <- array(0, dim = c(2,5))
  for (i in 1:2) {
    for (j in 1:5)  {
      logit <- glm(phenotypes[ ,i] ~ SNPs[ ,j], family=binomial(link=logit))
      odds_ratio <- exp(coef(logit))[2] # exponentiate coefficients to interpret as an odds ratio
      odds_ratios[i,j] <- odds_ratio
     }
  }
  return(odds_ratios) # return p_values in sorted order
}

"
> compute_odds()
          [,1]     [,2]      [,3]      [,4]      [,5]
[1,] 2.0746077 1.072721 0.7290407 0.8725775 0.9514117
[2,] 0.8550737 1.006670 1.6273190 0.8220741 0.9890246
"

# plots, part 5a
plot(phenotypes[ ,1], SNPs[ ,1]); abline(glm(phenotypes[ ,1] ~ SNPs[ ,1], family=binomial(link=logit)))
plot(phenotypes[ ,2], SNPs[ ,3]); abline(glm(phenotypes[ ,2] ~ SNPs[ ,3], family=binomial(link=logit)))


# Supporting function (For partial application of functions)
partial <- function(f, ...) {
  capture <- list(...)
  function(...) {
    do.call(f, c(capture, list(...)))
  }
}

# Colbert Seanker HW 5 Markov chain, Wright fisher model
minor.freq=400
alleles = 2*1000 # twice the number of individuals
Xn.1 <- function(Xn) rbinom(1,alleles, Xn/alleles)
chain <- function(X, count=1) {
    if (count==1500) return(X)
    c(X, 
      chain(Xn.1(X), count=count+1)
     )
}
png(file="Markov Chain trajectory", width=512, height=512,res=72)
plot(chain(minor.freq),xlab="time", ylab="number minor allele",
     main="Markov chain walk")
dev.off()
# The minor allele tends to vanish around 600-1000 steps
# The minor major allele rarely vanishes, but if it does, it does around 1000 steps

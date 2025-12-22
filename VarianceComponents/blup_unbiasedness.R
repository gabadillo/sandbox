# Chapter 3
set.seed(0)

# eq 39
s2a = 2
s2e = 2
n = 3
a = 1000

ai = matrix(rnorm(a, sd = sqrt(s2a)), ncol = 1)

BLUP = c()
BLUE = c()
r = 1000
for (i in 1:r){
  print(i)
  eij = matrix(rnorm(a*n, sd = sqrt(s2e)), nrow = a)
  yij = kronecker(matrix(1, nrow = 1, ncol = n), ai) + eij + 7
  
  yi. = matrix(apply(yij, 1, mean), ncol = 1)
  y.. = mean(yi.)
  
  ## ANOVA
  SSA = n*sum((yi. - y..)^2)
  MSA = SSA/(a-1)
  SSE = sum((yij - kronecker(matrix(1,nrow=1,ncol=n), yi.))^2)  
  MSE = SSE/(a*(n-1))
  
  # ANOVA Variance Components estimation 
  s2e = MSE
  s2a = (MSA - MSE)/n
  
  # BLUP
  ati = (MSA - MSE)/MSA*(yi. - y..)
  BLUP = c(BLUP, ati[1])
  BLUE = c(BLUE, yi.[1] - y..)
}

hist(BLUP); mean(BLUP); var(BLUP)
hist(BLUE); mean(BLUE); var(BLUE)


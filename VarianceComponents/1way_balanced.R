# Chapter 3
set.seed(0)

# eq 39
s2a = 1
s2e = 2
n = 10
a = 4

ai = matrix(rnorm(a, sd = sqrt(s2a)), ncol = 1)
eij = matrix(rnorm(a*n, sd = sqrt(s2e)), nrow = a)
yij = kronecker(matrix(1, nrow = 1, ncol = n), ai) + eij

yi. = matrix(apply(yij, 1, mean), ncol = 1)
y.. = mean(yi.)

## ANOVA
SSA = n*sum((yi. - y..)^2)
MSA = SSA/(a-1)
SSE = sum((yij - kronecker(matrix(1,nrow=1,ncol=n), yi.))^2)  
MSE = SSE/(a*(n-1))

# ANOVA hypothesis testing -> H0: s2a = 0
Ftest = MSA/MSE
pvalue = 1 - pf(Ftest, a-1, a*(n-1))

# ANOVA Variance Components estimation 
s2e = MSE
s2a = (MSA - MSE)/n

# BLUP
ati = (MSA - MSE)/MSA*(yi. - y..)


# lme4
library(lme4)
df = data.frame(y = c(yij), id = rep(1:a, n))
fm = lmer(y ~ (1|id), data = df)
uhat = fm@u
thetahat = as.data.frame(VarCorr(fm))[,'vcov']

pairs(cbind(ai, ati, uhat))

varcomps = rbind(c(s2a, s2e), thetahat)
colnames(varcomps) = c('s2a', 's2e')
rownames(varcomps) = c('ANOVA', 'lme4')
varcomps

# aov package
fm = aov(y ~ as.factor(id), data = df)
summary(fm)
Faov = summary(fm)[[1]]$`F value`[1]
# hypothesis testing
Faov; Ftest
# BLUEs 
aiaov = unname(fm$coefficients)
aiaov[-1] = aiaov[-1] + aiaov[1]
plot(aiaov, yi.); cor(aiaov, yi.)


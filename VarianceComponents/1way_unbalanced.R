# Chapter 3
set.seed(0)

# Section 3.6
s2a = 2
s2e = 2
a = 1000

ni = sample(1:3, a, replace = TRUE)
N = sum(ni)

ai = matrix(rnorm(a, sd = sqrt(s2a)), ncol = 1)

# matrix no longer valid, using lists instead
eij = map(ni, ~ rnorm(.x, mean = 0, sd = sqrt(s2e)))
yij = imap(eij, ~ .x + ai[.y])

yi. = sapply(yij, mean)
y.. = sum(ni*yi.)/N

## ANOVA
SSA = sum(ni*((yi. - y..)^2))
MSA = SSA/(a-1)
SSE = sum((unlist(yij) - rep(yi., ni))^2)
MSE = SSE/(N - a)

# ANOVA hypothesis testing -> H0: s2a = 0
Ftest = MSA/MSE
pvalue = 1 - pf(Ftest, a-1, N-a)

# ANOVA Variance Components estimation 
s2e = MSE
s2a = (MSA - MSE)/((N - sum(ni^2)/N)/(a-1))

# BLUP
hi = s2a/(s2a + s2e/ni)
ati = hi*(yi. - y..)


# lme4
library(lme4)
df = data.frame(y = unlist(yij), id = rep(1:a, ni))
fm = lmer(y ~ (1|id), data = df)
uhat = fm@u
thetahat = as.data.frame(VarCorr(fm))[,'vcov']

pairs(cbind(ai, ati, uhat))
plot(ati, uhat)

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


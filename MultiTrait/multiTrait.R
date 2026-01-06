library(BGLR)

# data
data(wheat)
X = wheat.X
K = AGHmatrix::Gmatrix(X)
Y = 100*wheat.Y
Y = cbind(Y, 2*Y[,1]+rnorm(nrow(Y), sd = 30))
rm(wheat.A, wheat.X, wheat.Y)


## Omega and R0
# UNS
CovUNS = list(type="UN")
# DIAG
CovDIAG = list(type="DIAG")
# FA
latent = 2
CovFA = list(type="FA", M = matrix(nrow = ncol(Y), ncol = latent, TRUE))

# Terms
B = list(X = X) # model = c('FIXED', 'BRR', 'SpikeSlab')
U = list(K = K, model = 'RKHS') # model = c('RKHS') 

# utils
n = nrow(X)
p = ncol(X)
t = ncol(Y)
J = matrix(1, nrow = n, ncol = 1)
H = matrix(1, nrow = t, ncol = 1)

##  Ex1 fit BRR with UNS and DIAG residual
A = B
A$model = 'BRR'
A$Cov = CovUNS
fm1 = Multitrait(y = Y, ETA = list(A = A), resCov = CovDIAG, 
                 nIter = 1200, burnIn = 200, saveAt = 'tmp/ex1_',
                 verbose = TRUE)

Omega = fm1$ETA$A$Cov$Omega
R0 = fm1$resCov$R
Beta = fm1$ETA$A$beta
yHat = fm1$ETAHat # = X%*%Beta + J%*%t(matrix(fm1$mu, ncol = 1))

##  Ex2 fit SpikeSlab with UNS and DIAG residual
A = B
A$model = 'SpikeSlab'
A$Cov = CovUNS
fm2 = Multitrait(y = Y, ETA = list(A = A), resCov = CovDIAG, 
                 nIter = 1200, burnIn = 200, saveAt = 'tmp/ex2_',
                 verbose = TRUE)

Omega = fm2$ETA$A$Cov$Omega
R0 = fm2$resCov$R
Beta = fm2$ETA$A$beta
yHat = fm2$ETAHat # = X%*%Beta + J%*%t(matrix(fm2$mu, ncol = 1))

##  Ex3 fit RKHS with FA(2) and DIAG residual
A = U
A$Cov = CovFA
fm3 = Multitrait(y = Y, ETA = list(A = A), resCov = CovDIAG, 
                 nIter = 1200, burnIn = 200, saveAt = 'tmp/ex3_',
                 verbose = TRUE)

W = fm3$ETA$A$Cov$W
PSI = fm3$ETA$A$Cov$PSI
Omega = fm3$ETA$A$Cov$Omega
R0 = fm3$resCov$R
u = fm3$ETA$A$u
yHat = fm3$ETAHat # = u + J%*%t(matrix(fm3$mu, ncol = 1))

# PCA vs AMMI vs FA

# PCA
Ys = scale(Y, scale = FALSE)
K = tcrossprod(t(Y))/(nrow(Y))

PCA = prcomp(K)
PCA = as.data.frame(PCA$x[,1:2])
PCA$EID = 1:5
PCA$Method = 'PCA'

# FA
FA = data.frame(W)
colnames(FA) = c('PC1', 'PC2')
FA$EID = 1:5
FA$Method = 'FA'

# AMMI
Yi. = matrix(apply(Y, 1, mean), ncol = 1)
Y.j = matrix(apply(Y, 2, mean), nrow = 1)
Y.. = matrix(mean(Y), ncol = 1, nrow = 1)


Z3 = Y - Yi.%*%t(H) - J%*%Y.j + J%*%Y..%*%t(H)
AMMI = svd(Z3)
AMMI = data.frame(AMMI$v[,1:2])
colnames(AMMI) = c('PC1', 'PC2')
AMMI$EID = 1:5
AMMI$Method = 'AMMI'

df = rbind(PCA, FA, AMMI)

library(tidyverse)
df %>% 
  ggplot(aes(x = PC1, y = PC2, color = as.factor(EID), label = EID))+
  geom_point(size = 5, alpha = 0.5)+
  facet_wrap(~Method, scales = 'free')
  


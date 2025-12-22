library(BGLR)
library(lme4GS)
library(sommer)
data(wheat)
set.seed(0)
# Gibbs sampler for mixed model y = µ + u + e 
# V(y) = V(a) + V(e) = Σ + R = Kσ2a + Iσ2e

# Gibbs sampler parameters
nIter = 4000
burnIn = 2000

# data
y = 2 * matrix(wheat.Y[,1], ncol = 1) + 100
K = tcrossprod(scale(wheat.X))/ncol(wheat.X)

## function
n = nrow(y)

mu.hist = matrix(nrow = 1, ncol = nIter)
u.hist = matrix(nrow = n, ncol = nIter)
varU.hist = matrix(nrow = 1, ncol = nIter)
varE.hist = matrix(nrow = 1, ncol = nIter)

R2 = 0.5
Sy = drop(var(y))
S0U = Sy*R2
S0E = Sy*(1 - R2)
df0 = 5

mu0 = 0
varMu = 1e10

mu = mean(y)
u = matrix(0, nrow = n, ncol = 1)
e = y - mu - u
varE = S0E
varU = S0U

I = diag(n)
Ki = solve(K + diag(n) * 1e-6)
df = n + df0 

for (i in 1:nIter){
  if (!i %% 50){print(i)}
  mu.hist[i] = mu
  u.hist[,i] = u
  varU.hist[i] = varU
  varE.hist[i] = varE
  
  # update µ
  e = e + mu
  C = n / varE + 1 / varMu
  RHS  = sum(e) / varE + mu0 / varMu
  mu = rnorm(1, RHS/C, sqrt(1/C))
  e = e - mu
  
  # update u
  e = e + u
  C = I / varE + Ki / varU
  RHS  = e / varE
  L = chol(C)
  w = backsolve(L, RHS, transpose = TRUE)
  uhat = backsolve(L, w, transpose = FALSE)
  delta = backsolve(L, rnorm(n), transpose = FALSE)
  u = uhat + delta
  e = e - u
  
  # update varU
  SS = t(u)%*%Ki%*%u + df0*S0U
  varU = c(SS/rchisq(1, df = df))
  
  # update varE
  SS = t(e)%*%e + df0*S0E
  varE = c(SS/rchisq(1, df = df))
}

plot(t(varE.hist))
plot(t(varU.hist))
plot(t(mu.hist))

muHat =  mean(mu.hist[-c(1:burnIn)])
varUhat = mean(varU.hist[-c(1:burnIn)])
varEhat = mean(varE.hist[-c(1:burnIn)])
uHat = apply(u.hist[,1:burnIn], 1, mean)
rHat = drop(cor(uHat, y))

mine_out = c('mu' = muHat, 'r' = rHat, 'varU' = varUhat, 'varE' = varEhat)

# BGLR
fm = BGLR(y = y, ETA = list(u = list(K = K, model = 'RKHS')), 
          nIter = nIter, burnIn = burnIn, saveAt = 'tmp/')
bglr_out = c('mu' = fm$mu, 'r' = cor(fm$ETA$u$u, y), 'varU' = fm$ETA$u$varU, 'varE' = fm$varE)
  
# LME4GS
ids = paste0('G', 1:n)
rownames(K) = colnames(K) = ids
df = data.frame(GID = ids, y = y)
out = lmerUvcov(y ~ (1|GID), data = df, Uvcov=list(GID = list(K=K)), verbose = TRUE)
# get random effects
tmp = ranef(out)[[1]]
uu = tmp[,1]
names(uu) = rownames(tmp)
uu = uu[df$GID]
# get varcomps
vv = as.data.frame(VarCorr(out))[,4]
lme4_out = c('mu' = unname(fixef(out)), 'r' = cor(uu, y), 'varU' = vv[1], 'varE' = vv[2])

# sommer
ans1 <- mmes(y~1,
             random= ~ vsm(ism(GID),Gu=K),
             rcov= ~ units, nIters=20,
             data=df, verbose = TRUE, 
             dateWarning = FALSE)
# get random effects
tmp = ans1$u
uu = tmp[,1]
uu = uu[df$GID]

mmer_out = c('mu' = unname(ans1$b[,1]), 'r' = cor(uu, y), 'varU' = ans1$theta[[1]], 'varE' = ans1$theta[[2]])


rbind(mine_out, bglr_out, lme4_out, mmer_out)

# mu         r     varU     varE
# mine_out 100.00007 0.8196479 2.137357 2.136970
# bglr_out  99.99899 0.8148897 2.136817 2.143441
# lme4_out 100.00000 0.8156659 2.118557 2.127987
# mmer_out 100.00000 0.8156643 2.118522 2.127999



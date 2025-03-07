#######################################################################################
#   Filename    :	DA_causal.R
#       
#   Required package  :    nleqslv, MASS, rjags, coda, ggplot2
#
########################################################################################

library("nleqslv")
library("MASS")
library("rjags")
library("coda")
library("ggplot2")
expit = function(x) exp(x)/(1+exp(x))

#####################################################
# Data preparation:
# Plug in the values for the variables defined below. 
# Note:
# Study s=1,...Z are with aggregated data only. 
# Study s=Z+1,...,K are with individual-level data
#####################################################

# n.cov: the number of covariates.
# n: the total number of patients across all studies.
# K: the total number of studies.
# Z: the total number of studies with aggregated data.

# mu.ag: a list recording the means of L in all K trials:
# mu.ag[[s]]: a numeric vector of length n.cov recording the mean of l1, l2, ... in study s. 

# v.ag: a list recording the sample variance of L in all K trials:
# v.ag[[s]]: a numeric vector of length n.cov recording the sample variance of l1, l2, ... in study s.

# y1.ag: a numeric vector of length K, where y1.ag[s] records E(Y|X=1,S=s).
# y1.n.ag: a numeric vector of length K, where y1.n.ag[s] records the number of patients with X=1,Y=1 in study s.
# y0.ag: a numeric vector of length K, where y0.ag[s] records E(Y|X=0,S=s).
# y0.n.ag: a numeric vector of length K, where y0.n.ag[s] records the number of patients with X=0,Y=1 in study s.

# y.ag: a numeric vector of length K, where y.ag[s] records the number of events in trial s.
# n.ag: a numeric vector of length K, where n.ag[s] records the sample size of trial s.

# n.xag: a numeric vector of length K, where n.xag[s] records the number of treated patients in trial s.

# ps: proportion of patients in each trial.
ps = n.ag/n
# px: randomization ratio in all trials, where px[s] records P(X=1|S=s).

# nj: a list of Z elements.
# nj[[z]] is a numeric vector of 2 elements, recording the number of treated and control patients in trial z, z=1...Z.

# mu: a list of Z elements.
# mu[[z]] is a matrix of dimension (n.cov + 1) * 2
# mu[[z]][,x+1] records the mean of L1, L2, ... and Y in treatment group x of trial z with aggregated data, where x=0,1 and z=1...Z.


# logrrjj: a numeric vector of K elements, where logrrjj[s] records the unstandardized log RR in study s.
# var.logrrjj: a numeric vector of K elements, where var.logrrjj[s] records the variance of the unstandardized log RR in study s.

# mu.q: a list of Z elements.
# mu.q[[z]] is a matrix of dimension (n.cov + 1) * 2.
# mu.q[[z]][,x+1] records the mean of L1^2, L2^2, ... and Y^2 in treatment group x of trial z with aggregated data, where x=0,1 and z=1...Z.

# data: a dataframe recording the IPD from trials s=Z+1, ..., K
#------------------------------------------------------------------------------------------------------
#            data = data.frame(int = 1, 
#                             x = df$'treatment', 
#                             l1 = df$'covariate 1', 
#                             l2 = df$'covariate 2',
#                             l3 = df$'covariate 3',
#                             ...
#                             ln = df$'covariate n',
#                             y = df$'outcome', 
#                             s = df$'study id') # Note: s can only take values Z+1, ..., K
#--------------------------------------------------------------------------------------------------------

# ---------------------------
# Step 1.1: Weight estimation
# ---------------------------

beta.func = function(est, AG, S, data){
  l.cov = paste0(rep("l",n.cov),1:n.cov)
  L = as.matrix(cbind(int = 1, data[,l.cov]))
  st = (data$s == S)*exp(L%*% est)
  indi = sapply(1:ncol(L), function(i) st*L[,i])
  part1 = colMeans(indi)*(dim(data)[1])/n
  part2 = c(1,mu.ag[[AG]])*(n.ag[AG]/n)
  return(part1-part2)
}

bstart = rep(0, n.cov + 1)
b = lapply((Z+1):K, function(k1)
  sapply(1:K, function(j1) nleqslv(bstart, beta.func, AG = j1, S = k1, data = data)$x)) 

# -------------------------------
# Step 1.2: effect standardization
# -------------------------------

# Truncation function
truncate = function(w_) {
  q3 = quantile(unlist(w_),.95)
  if (max(unlist(w_))>10) {
    a = w_
    a[which(a>=q3)] = q3
  } else {a = w_}
  return (as.vector(a))
}

ipw.m = function(j, k, data){
  l.cov = paste0(rep("l",n.cov),1:n.cov)
  L = as.matrix(cbind(int = 1, data[,l.cov]))
  data$w = as.vector(exp(L %*% b[[k-Z]][,j]))
  data$w = truncate(data$w)
  num = as.numeric(data$s == k) * data$y * data$x * data$w/(px[k]*ps[j])
  enum = as.numeric(data$s == k) * data$y * (1 - data$x) * data$w/((1-px[k])*ps[j])
  rrjk = (sum(num)/n) / (sum(enum)/n)
  return(rrjk)
}

rr = lapply((Z+1):K, function (ipd) sapply(1:K, function (ag) ipw.m(j = ag, k = ipd, data = data)))


# ---------------------------------
# Step 2.1: Pseudo-data generation
# ---------------------------------

l.cov = paste0(rep("l",n.cov),1:n.cov)
L = as.matrix(cbind(int = 1, data[,l.cov]))
# normal weights
w.ne.trunc = lapply((Z+1):K, function(k1) lapply(1:K, function(j1) as.vector(exp(L %*% b[[k1-Z]][,j1])))) 
# truncated weights
w = lapply(1:(K-Z), function(k1) lapply(1:K), function(j1) truncate(w0[[k1]][[j1]]))

# ----------- Step 2.1.1 Re-estimate E(L|S = j)------------------
mean.l = lapply(1:Z, function(j1) sapply((Z+1):K, function(k1) {
  try = (data$s == k1)*L*w[[k1-Z]][[j1]]
  return(colSums(try)/(n*ps[j1]))}))
mean.l = sapply(1:Z, function(j1) rowMeans(mean.l[[j1]]))

# ----------- Step 2.1.2 Estimate E(LaLb | S = j) -------------------

# A function to estimate E(u[,1]*u[,2] | S = j): 
cova.func = function(u){
  mu.cov = lapply(1:Z, function(j){
    sapply(1:(K-Z), function(k){
        sum((data$s == (k+Z)) * data[,u[1]] * data[,u[2]] * w[[k]][[j]])/(n*ps[j])
      })
    })
  mu.cov = lapply(1:Z, function(j) mean(mu.cov[[j]]))
  return(mu.cov)
}
ind2 = rep(1:n.cov,n.cov)
ind1 = ind2[order(ind2, decreasing = F)]
ind = data.frame(ind1 = ind1, ind2 = ind2)
ind$cov1 =  paste0("l",ind[,1])
ind$cov2 =  paste0("l",ind[,2])

# mu.lab$sj.x: E(LaLb|x,j) where La is ind[i,3] and Lb is ind[i,4]
mu.lab = lapply(1:dim(ind)[1], function(i) {
  u = c(ind[i,3],ind[i,4])
  cova.func(u = u)})
mu.lab = as.data.frame(do.call(rbind, mu.lab))
ag.ind = 1:Z
ag.ind = paste0(rep("s",Z),ag.ind)
colnames(mu.lab) = paste(ag.ind,sep=".")
mu.lab = cbind(ind, mu.lab)

# ----------- Step 2.1.3. Regenerate L1, L2 in trial j ----------
psd = lapply(1:Z, function(j){
    m.yl = mu[[j]][,1]*(1-px[j]) + mu[[j]][,2]*px[j]
    m.l = m.yl[1:n.cov]
    
    col.var = paste0("s",j)
    mat.l = matrix(as.numeric(mu.lab[,col.var]), nrow = n.cov, ncol = n.cov) 
    mat.l = mat.l - m.l %*% t(m.l)
    
    psd.j = as.data.frame(mvrnorm(n = sum(nj[[j]]),
                                  mu = m.l,
                                  Sigma = mat.l,
                                  empirical = TRUE)) # Pseudo-data of L
    
    psd.j$s = j
    psd.j$x = c(rep(0,nj[[j]][1]), rep(1,nj[[j]][2])) 
    psd.j$y = 0
    psd.j$y[which(psd.j$x == 0)][1:y0.n.ag[j]] = 1
    psd.j$y[which(psd.j$x == 1)][1:y1.n.ag[j]] = 1
    
    return(psd.j)})
psd = do.call(rbind, psd)
new = rbind(data, psd)

# -----------------------
# Step 2.2. M-estimation 
# -----------------------

# --------- Step 2.2.1. Estimating function -------

func = function(data, k1, ps, b.est, theta){
  
  # ps
  f0 = sapply(1:K, function(i) (data$s==i) - ps[i]) # P(S = j)
  #sapply(1:K, function(i) (data$s==i)*data$x/ps[i] - px[i])) # P(X=1|S) 
  
  # b.jk
  name.l = paste0(rep("l",n.cov), 1:n.cov)
  L = as.matrix(cbind(1, data[,name.l]))
  resi.s = lapply(1:K, function(j1)
    (data$s==k1)*exp(L %*% b.est[[k1-Z]][,j1])) # resi.s[[j]]: from k1 to source j1
  f2 = lapply(1:K, function(j1)
    t(sapply(1:ncol(L), function(i) resi.s[[j1]] * L[,i] - (data$s == j1)*L[,i]))) # f2[[j1]]: estimating function for b_j1k1 where k1 is fixed
  f2.out = do.call(rbind,f2)
  
  # theta
  f3 = sapply(1:K, function(j1)
    (data$s == k1) * data$y * exp(L %*% b.est[[k1-Z]][,j1]) * (data$x/px[j1] - theta[[k1-Z]][j1]*(1 - data$x)/(1-px[j1])))
  
  # summarize
  return(rbind(t(f0), f2.out, t(f3)))
}

# --------- Step 2.2.2. Matrix B ------------------

phi1 = lapply((Z+1):K, function(ipd) func(data = new, k1 = ipd, ps = ps, b.est = b, theta = rr))
mat.b1 = lapply((Z+1):K, function(ipd) var(t(phi1[[ipd-Z]])))

# --------- Step 2.2.3. Matrix A ------------------
# ps part
h = 1e-7
identity = diag(length(ps))
#mat.psx[[k]]: derivative wrt ps for the source trial k
mat.psx = lapply((Z+1):K, function(ipd){
  dev = lapply(1:length(ps), function(u){
    p. = complex(real = ps, imaginary = h*identity[u,])
    out = rowSums(Im(func(data = new,
                          k1 = ipd,
                          ps=p., 
                          b.est=b, 
                          theta=rr))/h)
    return(out)})
  mat.a.psx = -do.call(cbind,dev)/n
  return(mat.a.psx)
})

# b.jk
identity = diag(length(b[[1]][,1]))
dev = lapply((Z+1):K, function(ipd){
  lapply(1:K, function(j1){
    sapply(1:length(b[[1]][,1]), function(i){
      b.x= b
      b.x[[ipd - Z]][,j1] = complex(real = b[[ipd - Z]][,j1], imaginary = h*identity[i,])
      out = rowSums(Im(func(data = new, 
                            k1 = ipd,
                            ps=ps,
                            b.est=b.x, 
                            theta=rr))/h)
      return(out)})
  })
})
dev = lapply(1:(K-Z), function(i) do.call(cbind,dev[[i]]))
for (ipd in 1:(K-Z)) dev[[ipd]] = dev[[ipd]]/n
mat.bx = dev

# theta
identity = diag(length(rr[[1]]))
dev = lapply((Z+1):K, function(ipd){
  sapply(1:length(rr[[1]]), function(i){
    rr.x= rr
    rr.x[[ipd - Z]] = complex(real = rr.x[[ipd - Z]], imaginary = h*identity[i,])
    out = rowSums(Im(func(data = new, 
                          k1 = ipd,
                          ps=ps,
                          b.est=b, 
                          theta=rr.x))/h)
    return(out)})
})

for (ipd in 1:(K-Z)) dev[[ipd]] = dev[[ipd]]/n
mat.thetax = dev

mat.a = lapply(1:(K-Z), function(k1) cbind(mat.psx[[k1]], mat.bx[[k1]], mat.thetax[[k1]])) 

# ---------- Step 2.2.4. Variance estimates -----------

var.para1 = lapply(1:(K-Z), function(k1) (solve(mat.a[[k1]]) %*% mat.b1[[k1]] %*% t(solve(mat.a[[k1]])))/n)
n0.te = K + K*K + 1
n1.te = dim(var.para1[[1]])[1]

var.te = lapply(1:(K-Z), function(k1) var.para1[[k1]][n0.te:n1.te, n0.te:n1.te])

# ---------- Step 2.2.5. Delta method for log RR --------------

logrr = lapply(1:(K-Z), function(k1) log(rr[[k1]]))

var.logrr = lapply(1:(K-Z), function(k1) {
  delta = diag(1/rr[[k1]])
  return(delta %*% var.te[[k1]] %*% delta)
})

# ---------- Step 2.2.5. 95% CI of log RR ---------------

#ci[[k]][j,]: 95%CI of rr(j,k)
ci = lapply(1:(K-Z), function(k1) 
  t(sapply(1:K, function(j1) logrr[[k1]][j1] + qnorm(c(.025, .975))*sqrt(diag(var.logrr[[k1]])[j1]))))

# ----------------------------------------
# Step 3. Random-effect model estimation 
# ----------------------------------------

dat <- list("K" = K, "Z" = Z, 
            "y1" = t(do.call(rbind, logrr)), "v1" = do.call(cbind, var.logrr),
            "y2" = logrrjj, "v2" = diag(var.logrrjj))  

inits <- list(mu=0.0, 
              mu1=rep(0.0,K-Z), 
              s.tauj=1, 
              s.tau=1, 
              tau=rep(0,K),
              alpha=rep(0,K))

jagresult <- jags.model('.../metaRE14.txt',
                        data=dat, inits=inits, n.chains=5, n.adapt = 100000, quiet = T)
thesamps <- coda.samples(jagresult,c('mu','mu1','tau.sq','tauj.sq','tauk.sq','mu.pop'),n.iter=50000,thin=50)

# mu: summary treatment effect across populations
# mu.pop[j]: population j-specific treatment effect
# tau.sq: total heterogeneity
# tauj.sq: case-mix heterogeneity
# tauk.sq: residual heterogeneity (beyond case-mix heterogeneity and correlation)
final = summary(thesamps)$'quantiles'[1:(K+1),c(1,3,5)]

# ----------------------------
# Step 4. Drawing forest plot 
# ----------------------------
sample_data <- data.frame(study = c('Classic MA result', paste0("Population ",1:K)), 
                          index = 1:(K+1), 
                          result = final[,2], 
                          error_lower = final[,1], 
                          error_upper = final[,3]) 

te.rr = c(logrrjj[1:Z], unlist(logrr))
ci.jj = lapply(1:Z, function(k1) 
  logrrjj[k1] + qnorm(c(.025, .975))*sqrt(var.logrrjj[k1]))


ix.j = paste0(rep("j=", K*(K-Z)), c(1:Z, rep(1:K,K-Z)))
ix.k = paste0(rep(", k=", K*(K-Z)), c(1:Z, sort(rep((Z+1):K,K), decreasing = F)))

sample2 = data.frame(study = paste0(ix.j, ix.k),
                     index = 1:length(te.rr), 
                     logrr = te.rr, 
                     ci = rbind(do.call(rbind, ci.jj), do.call(rbind, ci)))   

# ------ Step 4.1. Standardized treatment effects --------
ggplot(data=sample2, aes(y=index, x=logrr, 
                         xmin=ci.1,  
                         xmax=ci.2)) + 
  geom_point() +  
  geom_errorbarh(height=.1) + 
  scale_y_continuous(breaks = 1:(K*(K-Z)+Z), labels=sample2$study) + 
  labs(title='Standardized treatment effects', 
       x='log relative risk and 95% confidence interval', 
       y = '')+ 
  geom_vline(xintercept=0, color='red', linetype='dashed', alpha=.8)

# ------ Step 4.2. Population-specific treatment effects ------------
ggplot(data=sample_data, aes(y=index, x=result, 
                             xmin=error_lower,  
                             xmax=error_upper)) + 
  geom_point() +  
  geom_errorbarh(height=.1) + 
  scale_y_continuous(labels=sample_data$study)+ 
  labs(title='Population-specific summary treatment effects', x='log relative risk and 95% credible interval', y = 'Population')+ 
  geom_vline(xintercept=0, color='red', linetype='dashed', alpha=.8)

# --------------- NUMERICAL RESULTS ----------------------

# random-effect model
print('Two-random-effect model')
(summary(thesamps))

# standardized treatment effects across studies
print('Standardized treatment effect estimates and 95%CI')
(sample2)

# classical meta-analysis results and population-specific summary treatment effects (estimate + 95%CI)
print('Population-specific summary treatment effects')
(sample_data)


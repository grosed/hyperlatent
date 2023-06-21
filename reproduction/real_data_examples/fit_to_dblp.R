library(praxi)

#source('../mcmc/getdat.R')
#source('../mcmc/mcmc.R')
#source("../mcmc/graphsummary.R") 

set.seed(890235)

####################### read in data #############################

load("./data/dblp_rw_to_fit.RData")
elst = rw_sub_elst
elst = lapply(elst, sort)
elst = unique(elst)
N = length(unique(unlist(elst)))
K = max( unlist( lapply( elst, length ) ) )

elen = unlist(lapply(elst, length))
dens = rep(NA, K-1)
for (k in 2:K){dens[k-1] = sum(elen==k)/choose(N,k) }

## ############################################
## 1) fit nonsimplicial model 
## ############################################

avec = rep(2, K-1)
bvec = rep(20, K-1)
lvec = .9*dens 
prior = list( m = c(0,0), sigmu = .2*diag(2), a0=avec, a1=avec, b0 = bvec, b1=bvec, lims=lvec, Phi = NA, nu = NA, lam = rep(1, K-1))
nIts = 400000
eps = c(0.15, .01) 

## run mcmc and save output
t1 = Sys.time()
out = runMCMC_cech(elst, N, K, nIts, eps, prior, modif=TRUE, bkids=c(39,71), DA=TRUE, rSimp=FALSE) 
t2 = Sys.time()
time = difftime(t2,t1, units="mins")

## save MCMC samples and timings 
save( out, file="./output/dblp_rw_nonsimp_fit.RData")
save( time, file="./output/dblp_rw_nonsimp_time.RData")

## ############################################
## 2) fit 'simplicial' model
## ############################################

avec = rep(2, K-1)
bvec = rep(20, K-1)
lvec = .9*dens
prior = list( m = c(0,0), sigmu = .2*diag(2), a0=avec, a1=avec, b0 = bvec, b1=bvec, lims=lvec, Phi = NA, nu = NA, lam = rep(1, K-1))
nIts = 400000
eps = c(0.15, .01) 

## run mcmc and save output
t1 = Sys.time()
out = runMCMC_cech(elst, N, K, nIts, eps, prior, modif=TRUE, bkids=c(39,71), DA=TRUE, rSimp=TRUE) 
t2 = Sys.time()
time = difftime(t2,t1, units="mins") 

## save MCMC samples and timings
save( out, file="./output/dblp_rw_simp_fit.RData")
save( time, file="./output/dblp_rw_simp_time.RData")

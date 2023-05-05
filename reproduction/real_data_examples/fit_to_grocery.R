source('../mcmc/mcmc.R')
source('../mcmc/getdat.R')
source("../mcmc/graphsummary.R") 

set.seed(1672348)

####################### read in data #############################

load("./data/gro_rndm_to_fit.RData") 
elst = rndm_sub_elst # to fit on
elst = lapply( elst, sort )
elst = unique( elst )
elen = unlist(lapply(elst, length))

## set up mcmc
N = max( unlist( elst ) )
K = max( unlist( lapply(elst, length) ) )
dens = rep(NA, K-1)
for( k in 2:K){ dens[k-1] = sum(elen==k)/choose(N, k) }

## ############################################
## 1) fit nonsimplicial model 
## ############################################

avec = rep(2,K-1)
bvec = rep(50,K-1)
limvec = .9*dens
prior = list( m = c(0,0), sigmu = .2*diag(2), a0=avec, a1=avec, b0 = bvec, b1=bvec, lims=limvec, Phi = NA, nu = NA, lam = rep(1, K-1))
nIts = 100000
eps = c(0.25, 0.25)

## run mcmc and save output
t1 = Sys.time()
out = runMCMC_cech(elst, N, K, nIts, eps, prior, bkids=c(1,9), modif=TRUE, DA=TRUE, rSimp=FALSE, rInc=TRUE)
t2 = Sys.time()
time = difftime(t2,t1,units="mins")

## save fitted model
save(out, file="./output/grocery_fit_nonsimp_rndm.RData")
save(time, file="./output/grocery_time_nonsim_rndm.RData")

## ############################################
## 2) fit 'simplicial' model
## ############################################

avec = rep(2,K-1)
bvec = rep(50,K-1)
limvec = .9*dens
prior = list( m = c(0,0), sigmu = .2*diag(2), a0=avec, a1=avec, b0 = bvec, b1=bvec, lims=limvec, Phi = NA, nu = NA, lam = rep(1, K-1))
nIts = 100000 
eps = c(.25, .25)

## run mcmc and save output
t1 = Sys.time()
out = runMCMC_cech(elst, N, K, nIts, eps, prior, bkids=c(1,9), modif=TRUE, DA=TRUE, rSimp=TRUE)
t2 = Sys.time()
time = difftime(t2,t1,units="mins")

## save fitted model
save(out, file="./output/grocery_fit_simp_rndm.RData")
save(time, file="./output/grocery_time_sim_rndm.RData")

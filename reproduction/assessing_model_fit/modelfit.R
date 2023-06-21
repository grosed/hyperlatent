library(hyperlatent)
library(MASS)

set.seed(258352)
N = 50 
K = 3
r = c(.35, .45) 
phi0 = rep(.0001, 2)
phi1 = rep(.0001, 2)
mu = c(0,0)
sigma = diag(2)
H = genHypergraph( N, K, r, phi0, phi1, mu, sigma, add_noise = TRUE, add_noise_full=FALSE, bkstn=FALSE)
prior = list( m = c(0,0), sigmu = .2*diag(2), a0=rep(.1,2), a1=rep(.1,2), b0 = rep(9.9,2), b1=rep(9.9,2), Phi = NA, nu = NA, lam = rep(1, K-1), lims=rep(1, K-1))
nIts = 115000
eps = c(0.15, 0.1)

nBrn = 75000  
ids = seq(nBrn, nIts, by=20)
print(paste("nIts = ", nIts, ", nsmp = ", length(ids)))

################ run MCMC ################

set.seed(266)
t1 = Sys.time()
mcmc <- runMCMC_cech(H$elst, N, K, nIts, eps, prior, modif=FALSE, DA=TRUE)
t2 = Sys.time() 
time = difftime(t2, t1, units="mins")

## plot the output 
strg = mcmc$strg
save( strg, file="./output/postfit.RData")
save( time, file="./output/post_time.RData")
save( H, file="./output/H.RData")

load( file="./output/postfit.RData" )
load( file="./output/H.RData" )

################ calculate predictives ################

post_dd = list(old=list(), oldk2=list(), oldk3=list(), new=list(), newk2=list(), newk3=list(), btwn=list(), btwnk2=list(), btwnk3=list())
true_dd_smp = list(old=list(), oldk2=list(), oldk3=list(), new=list(), newk2=list(), newk3=list(), btwn=list(), btwnk2=list(), btwnk3=list())
prmsTru = list(N = N, K = K, us = H$us, mu = H$mu, sigma = H$sigma, r = r, phi0 = phi0, phi1 = phi1)

ids = seq(nBrn, nIts, by=20)
nrep = length(ids)

Nstr = 5 # additional nodes

for( rep in 1:nrep ){

    ## SAMPLE FROM POSTERIOR PREDICTIVE ##
    
    ## sample hypergraph, gives eOld, eNew, eBtwn
    prms = list(N = N, K = K, us = strg$uTrc[,,ids[rep]], mu = strg$muTrc[,ids[rep]], sigma = strg$sigTrc[,,ids[rep]], r = strg$rTrc[,ids[rep]], phi0 = strg$phi0Trc[,ids[rep]], phi1 = strg$phi1Trc[,ids[rep]])
    htmp = genHyperG_with_NewConnection(Nstr, prms, periph=FALSE, addnoise=TRUE) 

    ## calculate degree distribution for eOld, eNew, eBtwn
    elen = unlist(lapply(htmp$eOld, length))
    post_dd$old[[rep]] = getDD_from_elist(htmp$eOld, N)
    post_dd$oldk2[[rep]] = getDD_from_elist(htmp$eOld[elen==2], N)
    post_dd$oldk3[[rep]] = getDD_from_elist(htmp$eOld[elen==3], N)

    htmp$eNew = lapply(htmp$eNew, function(x){ x - N } ) # index from 1:Nstr
    elen = unlist(lapply(htmp$eNew, length))
    post_dd$new[[rep]] = getDD_from_elist(htmp$eNew, Nstr) 
    post_dd$newk2[[rep]] = getDD_from_elist(htmp$eNew[elen==2], Nstr)
    post_dd$newk3[[rep]] = getDD_from_elist(htmp$eNew[elen==3], Nstr)

    elen = unlist(lapply(htmp$eBtwn, length))
    post_dd$btwn[[rep]] = getDD_from_elist(htmp$eBtwn, N+Nstr)
    post_dd$btwnk2[[rep]] = getDD_from_elist(htmp$eBtwn[elen==2], N+Nstr)
    post_dd$btwnk3[[rep]] = getDD_from_elist(htmp$eBtwn[elen==3], N+Nstr)

    ## SAMPLE FROM TRUE PARAMETERS ##

    ## sample hypergraph, gives eOld, eNew, dBtwn
    htmp = genHyperG_with_NewConnection(Nstr, prmsTru, periph=FALSE, addnoise=TRUE) 

    ## calculate degree distribution for eOld, eNew, eBtwn
    elen = unlist(lapply(htmp$eOld, length))
    true_dd_smp$old[[rep]] = getDD_from_elist(htmp$eOld, N)
    true_dd_smp$oldk2[[rep]] = getDD_from_elist(htmp$eOld[elen==2], N)
    true_dd_smp$oldk3[[rep]] = getDD_from_elist(htmp$eOld[elen==3], N)

    htmp$eNew = lapply(htmp$eNew, function(x){ x - N } ) # index from 1:Nstr
    elen = unlist(lapply(htmp$eNew, length))
    true_dd_smp$new[[rep]] = getDD_from_elist(htmp$eNew, Nstr)
    true_dd_smp$newk2[[rep]] = getDD_from_elist(htmp$eNew[elen==2], Nstr)
    true_dd_smp$newk3[[rep]] = getDD_from_elist(htmp$eNew[elen==3], Nstr)

    elen = unlist(lapply(htmp$eBtwn, length))
    true_dd_smp$btwn[[rep]] = getDD_from_elist(htmp$eBtwn, N)
    true_dd_smp$btwnk2[[rep]] = getDD_from_elist(htmp$eBtwn[elen==2], N+Nstr)
    true_dd_smp$btwnk3[[rep]] = getDD_from_elist(htmp$eBtwn[elen==3], N+Nstr)
    
    ## save every 100 iterations
    if (rep %% 100 == 0){
        save(post_dd, file="./output/postout_dd.RData")
        save(true_dd_smp, file="./output/true_dd_smp.RData")
    }
}
## save after last iteration
save(post_dd, file="./output/postout_dd.RData")
save(true_dd_smp, file="./output/true_dd_smp.RData") 

post_dd_av = list()
## old 
post_dd_av$old = getavDD_from_ddlist(post_dd$old)
post_dd_av$oldk2 = getavDD_from_ddlist(post_dd$oldk2)
post_dd_av$oldk3 = getavDD_from_ddlist(post_dd$oldk3)
## new
post_dd_av$new = getavDD_from_ddlist(post_dd$new)
post_dd_av$newk2 = getavDD_from_ddlist(post_dd$newk2)
post_dd_av$newk3 = getavDD_from_ddlist(post_dd$newk3)
## between
post_dd_av$btwn = getavDD_from_ddlist(post_dd$btwn)
post_dd_av$btwnk2 = getavDD_from_ddlist(post_dd$btwnk2)
post_dd_av$btwnk3 = getavDD_from_ddlist(post_dd$btwnk3)

save(post_dd_av, file="./output/post_dd_av.RData")

true_dd_smp_av = list()
## old 
true_dd_smp_av$old = getavDD_from_ddlist(true_dd_smp$old)
true_dd_smp_av$oldk2 = getavDD_from_ddlist(true_dd_smp$oldk2)
true_dd_smp_av$oldk3 = getavDD_from_ddlist(true_dd_smp$oldk3)
## new
true_dd_smp_av$new = getavDD_from_ddlist(true_dd_smp$new)
true_dd_smp_av$newk2 = getavDD_from_ddlist(true_dd_smp$newk2)
true_dd_smp_av$newk3 = getavDD_from_ddlist(true_dd_smp$newk3)
## between
true_dd_smp_av$btwn = getavDD_from_ddlist(true_dd_smp$btwn)
true_dd_smp_av$btwnk2 = getavDD_from_ddlist(true_dd_smp$btwnk2)
true_dd_smp_av$btwnk3 = getavDD_from_ddlist(true_dd_smp$btwnk3)

save(true_dd_smp_av, file="./output/true_dd_smp_av.RData")

## calculate true degree distribution
true_dd = list()
elen = unlist(lapply(H$elst, length))
true_dd$full = getDD_from_elist(H$elst, N)
true_dd$fullk2 = getDD_from_elist(H$elst[elen==2], N)
true_dd$fullk3 = getDD_from_elist(H$elst[elen==3], N)

save(true_dd, file="./output/true_dd.RData")

library(praxi)

## have 3 different types of misspec:
## 1) latent coordinates (normal, normal mixture, poisson process)
## 2) error (homog, non homog)
## 3) radii (homog, non homog)

## fixed parameters:
N = 50 
K = 3
d = 2

## read in id for misspec case
## id = as.numeric(commandArgs(trailingOnly=TRUE)[1])
## id = as.numeric(Sys.getenv('SGE_TASK_ID')) # task id, between 1 and 8
id = 1

## misspec key
umiss = c(1,2,3,4) # ie "nomiss", "clust1", "clust2", "unif"
rmiss = c(1,2) # ie "nomiss", "simp"
phimiss = c(1) # ie "nomiss"

cases = do.call( expand.grid, list(u=umiss, r=rmiss, phi=phimiss) )
case = cases[id, ] # misspec value of u, r and phi

## now:
## 1) generate data
## 2) run the mcmc
## 3) simulate from the pri and post predictives
ucase = case$u
rcase = case$r
phicase = case$phi
prms = list()
prms$N = N
prms$K = K

## set latent coordprms
if (ucase==1){ # no misspecification
    prms$mu = c(0,0)
    prms$sigma = .1*diag(d)
} else if (ucase==2){ # clust1 -> clear clusters
    prms$mu = matrix(c(0,0,.9,0), nrow=2)
    prms$sigma = array(c(c(.1*diag(d)),c(.1*diag(d))), dim=c(d,d,2))
    prms$lvec = c(.5, .5)
} else if (ucase==3){ # clust2 -> fuzzy clusters
    prms$mu = matrix(c(0,0,.5,0), nrow=2)
    prms$sigma = array(c(c(.1*diag(d)),c(.1*diag(d))), dim=c(d,d,2))
    prms$lvec = c(.5, .5)
} else if (ucase==4){ # uniform
    ## nothing needed for this case
}

## set radii parameters
if (rcase==1){ # no misspecification
    prms$r = c(.17, .19)
} else if (rcase==2){ # simplicial hypergraph (r equal)
    prms$r = c(.17, .17)
} else if (rcase==3){ # nonhom mn
    ## mn is determined by scaling
    rvec = c(.17, .19)
    mnr = mean(rvec)
    diff2mn = mnr - rvec[1]
    prms$ravec = c(10,10)
    prms$rbvec = c(10,10)
    prms$rlim = c( rvec[1] - diff2mn, mnr, rvec[2] + diff2mn)
} else if (rcase==4){ # nonhom bndry
    rvec = c(.17, .19)
    mnr = mean(rvec)
    diff2mn = mnr - rvec[1]
    prms$ravec = c(.25,.25)
    prms$rbvec = c(.25,.25)
    prms$rlim = c( rvec[1] - diff2mn, mnr, rvec[2] + diff2mn)
}

## set noise parameters
if (phicase==1){ # no misspecificaion
    prms$phi0 = c(.0001, .0001)
    prms$phi1 = c(.0001, .0001)
} else if (phicase==2){ # phi from mean
    prms$a0vec = c(5, 5)
    prms$b0vec = c(5, 5)
    prms$a1vec = c(5, 5)
    prms$b1vec = c(5, 5)
} else if (phicase==3){ # phi from boundary
    prms$a0vec = c(.5, .5)
    prms$b0vec = c(.5, .5)
    prms$a1vec = c(.5, .5)
    prms$b1vec = c(.5, .5)
}

####################### SIMULATE H #####################################

## set seed
set.seed(id+100)

## save parameter list
label = paste("prms_u",ucase,"r",rcase,"phi",phicase,sep="")
prms$label = label

## gen graph and save
ucase = c("nomiss", "clust", "clust", "unif")[ucase]
rcase = c("nomiss", "simp", "nonhom", "nonhom")[rcase]
phicase = c("nomiss", "nonhom", "nonhom")[phicase]

keep = FALSE 
while( keep == FALSE ){
    t1 = Sys.time()
    H = genHyp_shell( N, K, d, prms, ucase, rcase, phicase )
    t2 = Sys.time()

    ## make sure hypergraph is connected with N nodes
    H = get_largest_conn_comp(H, N)
    Ntmp = dim(H$us)[1]
    print( Ntmp )
    is_conn = check_connected(H$elst, Ntmp) 
    if( abs( diff( c(dim(H$us)[1], N) ) ) <= 0 & is_conn==TRUE ){keep = TRUE}
}

## save parameters and simulated hypergraphs
save( prms, file=paste("./output/", label,".RData", sep=""))
save( H, file=paste("./output/H_", label, ".RData",sep=""))

####################### RUN MCMC ########################################
elen = unlist(lapply(H$elst, length))
hdens = rep(0, K-1)
for( k in 2:K ){ hdens[k-1] =  sum(elen==k) / choose(N,k) }
prior = list( m = c(0,0), sigmu = .2*diag(2), a0=rep(.1,2), a1=rep(.1,2), b0 = rep(9.9,2), b1=rep(9.9,2), Phi = NA, nu = NA, lam = rep(1, K-1), lims=hdens)
nIts = 70000
nBrn = 35000
eps = c(.01, .01)
            
## fit the model for this data with DA since faster
t1 = Sys.time()
mcmcout <- runMCMC_cech(H$elst, N, K, nIts, eps, prior, DA=TRUE)
t2 = Sys.time()
timetot = difftime(t2, t1, units="hours")

strg = mcmcout$strg
## save mcmc output and timings
save( mcmcout, file=paste("./output/mcmcout_", label, ".RData", sep=""))
save( timetot, file=paste("./output/time_", label, ".RData", sep=""))

####################### SUMMARY OF DATA ###############################
## this is to compare with predictives

elen = unlist(lapply( H$elst, length))

truesum = list()
truesum$DD = getDD_from_elist(H$elst, N)
truesum$ntri = get_ntri( H$elst[elen==2] )
hypcounts = get_hypmotif_count( H$elst )
truesum$h1 = hypcounts[1]
truesum$h2 = hypcounts[2]
truesum$h3 = hypcounts[3]
truesum$ijk = sum( elen==3 )
truesum$nclst = get_nclst( H$us )

save( truesum, file=paste("./output/true_summary_", label, ".RData", sep=""))

####################### PREDICTIVE SUMMARY #######################################
id_pst = seq(nBrn, nIts, 10) 

## since these are model parameters, form does not change
pstTrc = list(uTrc = mcmcout$strg$uTrc[,,id_pst], muTrc = mcmcout$strg$muTrc[,id_pst], sigTrc = mcmcout$strg$sigTrc[,,id_pst], rTrc = mcmcout$strg$rTrc[,id_pst], phi0Trc = mcmcout$strg$phi0Trc[,id_pst], phi1Trc = mcmcout$strg$phi1Trc[,id_pst] ) 

## sample from predictive (what to record?)
preds = sim_from_preds(pstTrc, ucase, rcase, phicase, H, N, K, d=2)
save( preds, file=paste("./output/pred_summary_", label, ".RData", sep=""))

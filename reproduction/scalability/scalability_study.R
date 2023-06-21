library(hyperlatent)

set.seed(52365)

Nvec = c(10, 30, 50, 100) # this is approximate list of Ns (accept largest component close to)
nNs = length(Nvec)
K=3
d=2
mu=c(0,0)
phi = c(.001, .001)
r_sml = c(.75, .85)
r_lrg = c(1.4, 1.5)
svec = c(2, 2.5, 3, 4.5) # scale s, keep radii the same 

nRep = 10
dens = array(NA, c(2, nRep, K-1, nNs)) # r smaller/larger, reps, density, nodes
nmat = array(NA, c(2, nRep, nNs)) # rsmaller/larger, reps, nvalues
timings = array(NA, c(2, nRep, 2, nNs)) # r smaller/large, reps, DA/non, nodes

nIts = 200
eps_sml = c(0.1, .05)
eps_lrg = c(0.05, .05)

for (n in 1:nNs){

    for (rep in 1:nRep){

        ## r SMALLER ##
        
        ## sample connected hypergraph
        keep = FALSE
        while( keep == FALSE ){
            H = genHypergraph( Nvec[n], K, r_sml, phi, phi, mu, svec[n]*diag(d), add_noise = FALSE, bkstn=FALSE)
            H = get_largest_conn_comp(H, Nvec[n]) 
            nmat[1,rep,n] = dim(H$us)[1]
            is_conn = check_connected(H$elst, nmat[1,rep,n]) # should be connected
            if( abs( diff( c(dim(H$us)[1], Nvec[n]) ) ) <= 3 & is_conn==TRUE ){keep = TRUE}
        }
        print( paste( "Sampled H with r_sml, N = ", nmat[1,rep,n], " and rep = ", rep ) )
        
        ## record the density
        elen = unlist(lapply(H$elst, length))
        dens[1,rep,,n] = c( sum(elen==2)/choose(nmat[1,rep,n],2), sum(elen==3)/choose(nmat[1,rep,n],3) )

        ## run MCMC with DA
        prior = list( m = c(0,0), sigmu = .2*diag(2), a0=rep(1,K-1), a1=rep(1,K-1), b0 = rep(9,K-1), b1=rep(9,K-1), Phi = NA, nu = NA, lam = rep(1, K-1), lims=rep(1,K-1))
        tst = runMCMC_cech(H$elst, nmat[1,rep,n], K, nIts, eps_sml, prior, modif=FALSE, DA=TRUE)
        timings[1,rep,1,n] = tst$time / nIts

        ## run MCMC without DA
        prior = list( m = c(0,0), sigmu = .2*diag(2), a0=rep(1,K-1), a1=rep(1,K-1), b0 = rep(9,K-1), b1=rep(9,K-1), Phi = NA, nu = NA, lam = rep(1, K-1), lims=rep(1,K-1))
        tst = runMCMC_cech(H$elst, nmat[1,rep,n], K, nIts, eps_sml, prior, modif=FALSE, DA=FALSE)
        timings[1,rep,2,n] = tst$time / nIts

    }
    ## save current values given stage of smaller r study is completed
    save(dens, file="./output/scalability_dens.RData")
    save(nmat, file="./output/scalability_nmat.RData")
    save(timings, file="./output/scalability_times.RData")
    print( paste( "Finished smaller r for N roughly = ", Nvec[n] ) )

    for (rep in 1:nRep){
        ## r LARGER ##

        ## sample connected hypergraph
        keep = FALSE
        while( keep == FALSE ){
            H = genHypergraph( Nvec[n], K, r_lrg, phi, phi, mu, svec[n]*diag(d), add_noise = FALSE, bkstn=FALSE)
            H = get_largest_conn_comp(H, Nvec[n])
            nmat[2,rep,n] = dim(H$us)[1]
            is_conn = check_connected(H$elst, nmat[2,rep,n]) # should be connected
            if( abs( diff( c(dim(H$us)[1], Nvec[n]) ) ) <= 3 & is_conn==TRUE ){keep = TRUE}
        }
        print( paste( "Sampled H with r_lrg, N = ", nmat[2,rep,n], " and rep = ", rep ) )
        
        ## record the density
        elen = unlist(lapply(H$elst, length))
        dens[2,rep,,n] = c( sum(elen==2)/choose(nmat[2,rep,n],2), sum(elen==3)/choose(nmat[2,rep,n],3) )

        ## run MCMC with DA
        prior = list( m = c(0,0), sigmu = .2*diag(2), a0=rep(1,K-1), a1=rep(1,K-1), b0 = rep(9,K-1), b1=rep(9,K-1), Phi = NA, nu = NA, lam = rep(1, K-1), lims=rep(1,K-1))
        tst = runMCMC_cech(H$elst, nmat[2,rep,n], K, nIts, eps_lrg, prior, modif=FALSE, DA=TRUE)
        timings[2,rep,1,n] = tst$time / nIts

        ## run MCMC without DA
        prior = list( m = c(0,0), sigmu = .2*diag(2), a0=rep(1,K-1), a1=rep(1,K-1), b0 = rep(9,K-1), b1=rep(9,K-1), Phi = NA, nu = NA, lam = rep(1, K-1), lims=rep(1,K-1))
        tst = runMCMC_cech(H$elst, nmat[2,rep,n], K, nIts, eps_lrg, prior, modif=FALSE, DA=FALSE)
        timings[2,rep,2,n] = tst$time / nIts

    }
    ## save current values given stage of larger r study is completed
    save(dens, file="./output/scalability_dens.RData")
    save(nmat, file="./output/scalability_nmat.RData")
    save(timings, file="./output/scalability_times.RData")
    print( paste( "Finished larger r for N roughly = ", Nvec[n] ) )
}

## save everything 
save(dens, file="./output/scalability_dens.RData")
save(nmat, file="./output/scalability_nmat.RData")
save(timings, file="./output/scalability_times.RData")

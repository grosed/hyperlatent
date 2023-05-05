## script containing MCMC function, and plotting
source("../mcmc/mcmc_functions.R")
library('MCMCpack')
library('mvtnorm')
library('extraDistr')
library('ramcmc')

initStrg <- function( N, nIts, K, d )
{
                                        # function to intiialise storage
                                        # K is max e order, for {ij} graph K=2
    
    strg = list( uTrc = array( data=NA, dim=c(N,d,nIts) ), # stores coordinate trace
            muTrc = array( data=NA, dim=c(d, nIts) ), # stores mu trace
            sigTrc = array( data=NA, dim=c(d,d, nIts ) ), # stores sigmamu trace
            phi0Trc = array( data=NA, dim=c(K-1,nIts) ), # stores hyperedge noise
            phi1Trc = array( data=NA, dim=c(K-1,nIts) ), # stores hyperedge noise
            rTrc = array( data=NA, dim=c(K-1,nIts) ), # stores the radii
            ll = array( data=NA, dim=c(K-1,nIts) ) # stored log likelihood (one for each k)
            ) 

    return( strg )
}

####################### MCMC CODE ##########################

runMCMC_cech <- function(elstH, N, K, nIts, eps, prior, modif=TRUE, ulam=.8, bkids=NA, DA=TRUE, rSimp=FALSE, rInc=TRUE)
{
                                        # elstH is the list of observed edges
                                        # N is number of nodes in graph
                                        # K is maximum hyperedge order = length(r)+1
                                        # nIts is number of iterations
                                        # eps is vector of epsilon for rw MH eps[1] for u
                                        # prior contains prior info
                                        # modif == TRUE indicates to use modified model
                                        # DA == TRUE for delayed acceptance update
                                        # rSimp == TRUE if 'simplicial' variant so rk = rk-1
                                        # rInc == TRUE is rk > rk-1, otherwise ! rk < rk-1

    d = length( prior$m )
    ## Set up storage
    strg = initStrg( N, nIts, K, d )
    llPrp = rep(NA, K-1) # stores proposal log-likelihood
    count_u = rep(0,N-d) # acceptance rates for coordinates
    count_r = rep(0, K-1) # acceptance rates for radii
    S_u = eps[1]*array( diag(d), dim=c(d,d,N-d) )
    S_r = eps[2]*array( diag(1), dim=c(1,1,K-1) )
    da = rep(0, K-1) # monitor rejections in DA

    ## make sure observed hyperedges are ordered
    elstH = lapply( elstH, sort ) 
    
    ## INITIALISE ##

    ## Init latent coordinates
    ## get coords from GMDS
    uinit <- get_u_init(elstH, N, K, d, ulam)
    ## apply bookstein
    bk = getBookstein( uinit$coords, d, PCA=TRUE, bkids)
    fixed = bk$id # index of fixed coordinates
    print( paste( "bookstein coords = ", fixed))
    strg$uTrc[,,1] = bk$u

    ## Init phi0 and phi1
    if (modif == TRUE ){ #model with two noise prms
        for (k in 2:K){
            strg$phi0Trc[k-1,1] = rnsbeta(1, prior$a0[k-1], prior$b0[k-1], min=0, max=prior$lims[k-1] )
            strg$phi1Trc[k-1,1] = rnsbeta(1, prior$a1[k-1], prior$b1[k-1], min=0, max=prior$lims[k-1] )
        }
    } else { # both noise parameters are the same
        for (k in 2:K){
            strg$phi0Trc[k-1,1] = rnsbeta(1, prior$a0[k-1], prior$b0[k-1], min=0, max=prior$lims[k-1] )
            strg$phi1Trc[k-1,1] = strg$phi0Trc[k-1,1]
        }
    }
   
    ## Init radii
    strg$rTrc[,1] = get_r_init(strg$uTrc[,,1], elstH, N, K, rInc) 
    if (rSimp==TRUE){ strg$rTrc[,1] = rep( mean(strg$rTrc[,1], K-1 ) ) }    

    ## Init mu and sigma (and their priors)

    ## get sigma and mu prior
    prior$nu = 20
    prior$Phi = prior$nu * get_sig_prior( bk$u )
    prior$m = rep(0,d) 
    prior$sigmu = .1*diag(d)

    ## init sigma, mu 
    musiginit = get_mu_and_sig_init(elstH, N, K, prior, strg$rTrc[,1], strg$phi0Trc[,1], strg$phi1Trc[,1], fixed, d, tol=.5*choose(N,K))
    strg$muTrc[,1] = apply( musiginit$mus, 2, mean )
    strg$sigTrc[,,1] = apply( musiginit$covs, c(2,3), mean )
    if (0){ # if too slow, can init as:
        strg$muTrc[,1] = c(0,0)
        strg$sigTrc[,,1] = diag(d)
    }
    
    ## Induce graph and initalise likelihood
    if (rSimp==FALSE){
        elstCur = get_graph_from_cech(strg$uTrc[,,1], strg$rTrc[,1])
    } else if (rSimp==TRUE){
        elstCur = get_simplicial_graph_from_cech(strg$uTrc[,,1], strg$rTrc[1,1], K)
    }
    strg$ll[,1] = loglikelihood(elstH, elstCur, N, K, strg$phi0Trc[,1], strg$phi1Trc[,1] )

    ## RUN MCMC ##
    t1 = Sys.time()
    for (i in 2:nIts)
    {

        ## UPDATE LATENT COORDINATES U via MH ##
        strg$uTrc[,,i] = strg$uTrc[,,i-1] # current latent coordinates
        strg$ll[,i] = strg$ll[,i-1] # current log likelihood 

        if (DA == TRUE){ ## update latent coordinates with DA
            
            for (n in (1:N)[-fixed]) # loop over nodes(-fixed BK coords)
            {
                nd_id = which((1:N)[-fixed]==n) # find node index in looping vector
                
                ## propose new coordinate (random walk on nth node)
                uPrp = strg$uTrc[,,i]
		p_type = sample(1:4, 1, prob=c(.8, 0.2, 0, 0))
		if (p_type==1){
		   etmp = rnorm(d, mean=0, sd = 1) 
		   uPrp[n,] = uPrp[n,] + S_u[,,nd_id] %*% etmp
		} else if (p_type==2){
		   uPrp[n,] = uPrp[n,] + rnorm(d, mean=0, sd = 0.1)
		} else if (p_type==3){
		   uPrp[n,] = uPrp[n,] + rnorm(d, mean=0, sd = 0.01)
   		} else if (p_type==4){	  
		   uPrp[n,] = uPrp[n,] + rnorm(d, mean=0, sd = 0.001)
		}

                ## delayed acceptance 
		uprior = dmvnorm( matrix(uPrp[n,], nrow=1), mean=strg$muTrc[,i-1], sigma=strg$sigTrc[,,i-1], log=TRUE ) - dmvnorm( matrix(strg$uTrc[n,,i], nrow=1), strg$muTrc[,i-1], strg$sigTrc[,,i-1], log=TRUE ) # prior
		ARapp = 0
                for (k in 2:K){

                    ## calculate order k hyperedges graph induced from proposed coords
                    elstPrp = get_graph_from_cech_order_k(uPrp, strg$rTrc[,i-1], k)

                    ## calculate order k log likelihood
                    llPrp[k-1] = ll_order_k(elstH, elstPrp, N, k, strg$phi0Trc[,i-1], strg$phi1Trc[,i-1]) 

                    ## order k A/R ratio
		    if (k==2){
		       AR = uprior + llPrp[k-1] - strg$ll[k-1,i]
		    } else {
		       AR = llPrp[k-1] - strg$ll[k-1,i]	
		    }
		    ARapp = ARapp + AR
                    
                    ## early reject or continue?
                    if ( log(runif(1)) < AR ){
                        id_accpt = TRUE # continue to next DA
                    } else {
                        id_accpt=FALSE # reject and break
                        da[k-1] = da[k-1] + 1 # update counter for delayed acceptance
                        break
                    }
                }

                ## if accept proposal, update current state
                if (id_accpt == TRUE){
                    ## update 'current' state
                    strg$ll[,i] = llPrp
                    strg$uTrc[,,i] = uPrp       
                    ## update count
                    count_u[nd_id] = count_u[nd_id] + 1            
                }

		## update S_u if in burn-in
		if ( (i < (nIts/2)) & p_type==1 ){ S_u[,,nd_id] = adapt_S(S_u[,,nd_id], etmp, min(c(1,exp(ARapp))), i) }
            }

        } else if (DA == FALSE){ ## Update latent coordinates without DA

            for (n in (1:N)[-fixed]) # loop over nodes (-fixed BK coords)
            {
                nd_id = which((1:N)[-fixed]==n) # find node index in looping vector
                
                ## propose new coordinate (random walk on nth node)
                uPrp = strg$uTrc[,,i]
		p_type = sample(1:4, 1, prob=c(0.8, 0.2, 0, 0))
		if (p_type==1){
		   etmp = rnorm(d, mean=0, sd = 1) # sample from N(0, diag(2))
		   uPrp[n,] = uPrp[n,] + S_u[,,nd_id] %*% etmp
		} else if (p_type==2){
		   uPrp[n,] = uPrp[n,] + rnorm(d, mean=0, sd = 0.1)
		} else if (p_type==3){
		   uPrp[n,] = uPrp[n,] + rnorm(d, mean=0, sd = 0.01)
   		} else if (p_type==4){	  
		   uPrp[n,] = uPrp[n,] + rnorm(d, mean=0, sd = 0.001)
		}
                
                ## calculate graph induced from proposed coords (cheaper in simplicial case)
                if (rSimp==FALSE){
                    elstPrp = get_graph_from_cech(uPrp, strg$rTrc[,i-1])
                } else if (rSimp==TRUE){
                    elstPrp = get_simplicial_graph_from_cech(uPrp, strg$rTrc[1,i-1], K)
                }
                
                ## calculate log likelihood
                llPrp = loglikelihood(elstH, elstPrp, N, K, strg$phi0Trc[,i-1], strg$phi1Trc[,i-1])
                ## AR ratio
                dprior <- dmvnorm( uPrp[n,], strg$muTrc[,i-1], strg$sigTrc[,,i-1], log=TRUE ) - dmvnorm( strg$uTrc[n,,i], strg$muTrc[,i-1], strg$sigTrc[,,i-1], log=TRUE )
                AR <- sum(llPrp) - sum(strg$ll[,i]) + dprior

                if ( log(runif(1)) < AR )
                {
                    ## update 'current' state
                    strg$ll[,i] = llPrp 
                    strg$uTrc[,,i] = uPrp
                    
                    ## update count
                    nd_id = which((1:N)[-fixed]==n) # find node index in looping vector
                    count_u[nd_id] = count_u[nd_id] + 1            
                } 
            }
            
            ## update S_u if in burn-in
	    if ( (i < (nIts/2)) & p_type==1){ S_u[,,nd_id] = adapt_S(S_u[,,nd_id], etmp, min(c(1,exp(AR))), i) }
        }
        
        ## UPDATE RADII ##
        strg$rTrc[,i] = strg$rTrc[,i-1] # 'current' r value
     
        if (rSimp==FALSE){ # non-simplicial model with rk > rk-1
            for (k in 2:K){ 

                ## propose new value
                rPrp = strg$rTrc[,i]
		p_type = sample(1:3, 1, prob=c(1, 0, 0))
		if (p_type==1){
		   etmp = rnorm(1)
		   rPrp[k-1] = rPrp[k-1] + S_r[,,k-1] %*% etmp
		} else if (p_type==2){
		   rPrp[k-1] = rPrp[k-1] + rnorm(1, mean=0, sd = 0.1)
		} else if (p_type==3){
		   rPrp[k-1] = rPrp[k-1] + rnorm(1, mean=0, sd = 0.01)	
		}

                if (rInc==TRUE){
                    r_invld = any( diff(rPrp) <= 0 ) # reject r if TRUE
                } else if (rInc==FALSE){
                    r_invld = all( diff(rPrp) <= 0 ) # reject r if TRUE
                }

                if ( r_invld | sum( rPrp <=0 ) >= 1 ){
                    AR =  -1e6
                    print( "Reject proposal r" )
                } else {
                    ## calculate log likelihood 
                    elstPrp = get_graph_from_cech_order_k( strg$uTrc[,,i], rPrp, k)
                    llPrp[k-1] = ll_order_k(elstH, elstPrp, N, k, strg$phi0Trc[,i-1], strg$phi1Trc[,i-1])

                    ## calculate AR
                    dprior = dexp(rPrp[k-1], prior$lam[k-1], log=TRUE) - dexp(strg$rTrc[k-1,i], prior$lam[k-1], log=TRUE) 
                    AR = llPrp[k-1] - strg$ll[k-1,i] + dprior
                }
		
                ## accept/reject
                if( log(runif(1)) < AR )
                {
                    ## update 'current' state
                    strg$ll[k-1,i] = llPrp[k-1]
                    strg$rTrc[k-1,i] = rPrp[k-1]
                    
                    ## update count
                    count_r[k-1] = count_r[k-1] + 1
                }

                ## update S_r if in burn-in
		if ( (i < (nIts/2)) * (i > 100) & p_type==1 ){S_r[,,k-1] = adapt_S(S_r[,,k-1], etmp, min(c(1,exp(AR))), i)}
            }
        } else if (rSimp==TRUE){ # model with rk = rk-1

            ## propose new value
            rPrp = strg$rTrc[1,i]
	    p_type = sample(1:3, 1, prob=c(1, 0, 0))
	    if (p_type==1){
	       etmp = rnorm(1)
	       rPrp[1] = rPrp[1] + S_r[,,1] %*% etmp
	    } else if (p_type==2){
	       rPrp[1] = rPrp[1] + rnorm(1, mean=0, sd = 0.1)
	    } else if (p_type==3){
	       rPrp[1] = rPrp[1] + rnorm(1, mean=0, sd = 0.01)	
            }

            if ( rPrp <= 0 ){
                AR =  -1e6
                print( "Reject proposal r" )
            } else {

                ## calculate log likelihood
                elstPrp = get_simplicial_graph_from_cech(strg$uTrc[,,i], c(rPrp), K)
                for (k in 2:K){
                    llPrp[k-1] = ll_order_k(elstH, elstPrp, N, k, strg$phi0Trc[,i-1], strg$phi1Trc[,i-1])
                }
                ## calculate AR
                dprior = dexp(rPrp, prior$lam[1], log=TRUE) - dexp(strg$rTrc[1,i], prior$lam[1], log=TRUE) 
                AR = sum(llPrp) - sum(strg$ll[,i]) + dprior 
            }

            ## accept/reject
            if( log(runif(1)) < AR )
            {
                ## update 'current' state
                strg$ll[,i] = llPrp
                strg$rTrc[,i] = rep(rPrp, K-1)
                
                ## update count
                count_r[1] = count_r[1] + 1
            }
            ## update S_r if in burn-in
	    if ( (i < (nIts/2)) & (i > 100) & p_type==1 ){S_r[,,1] = adapt_S(S_r[,,1], etmp, min(c(1,exp(AR))), i)}
        }

        ## UPDATE NOISE PHI ##
        if (rSimp==FALSE){
            elstCur = get_graph_from_cech(strg$uTrc[,,i], strg$rTrc[,i])
        } else if (rSimp==TRUE){
            elstCur = get_simplicial_graph_from_cech(strg$uTrc[,,i], strg$rTrc[1,i], K)
        }
        for (k in 2:K)
        {
            ## calculate M11, M10, M01, M00 
            lenH = lengths( elstH )
            lenG = lengths( elstCur )
            nCmb = choose(N,k)
            Hktmp = elstH[ lenH == k ]
            Gktmp = elstCur[ lenG == k ]
	    M11 = dim(intersect( convert_l_to_df(Gktmp, k), convert_l_to_df(Hktmp, k)))[1]
            M10 = length( Gktmp ) - M11
            M01 = length( Hktmp ) - M11
            M00 = nCmb - (M11 + M10 + M01)

            ## sample \upvarphi_k^0 and \upvarphi_k^1
            if (modif==TRUE){
                strg$phi0Trc[k-1,i] = rnsbeta(1, prior$a0[k-1] + M01 , prior$b0[k-1] + M00, min=0, max=prior$lims[k-1] )
                strg$phi1Trc[k-1,i] = rnsbeta(1, prior$a1[k-1] + M10, prior$b1[k-1] + M11, min=0, max=prior$lims[k-1] )
            } else {
                strg$phi0Trc[k-1,i] = rnsbeta(1, prior$a0[k-1] + M01 + M10, prior$b0[k-1] + M00 + M11, min=0, max=prior$lims[k-1])
                strg$phi1Trc[k-1,i] = strg$phi0Trc[k-1,i]
            }
        }
        
        ## UPDATE MU via Gibbs ##
        invsig = solve(strg$sigTrc[,,i-1])
        invsigmu = solve(prior$sigmu)
        sig = solve( N * invsig + invsigmu )
        mn = sig %*% ( rowSums( invsig %*% t(strg$uTrc[,,i])) + invsigmu %*% prior$m )
        strg$muTrc[,i] = mvrnorm( n=1, mu=mn, Sigma=sig )

        ## UPDATE SIGMA via Gibbs ##
        diff = t( t(strg$uTrc[,,i]) - strg$muTrc[,i] ) 
        Phi = prior$Phi + t(diff) %*% diff
        strg$sigTrc[,,i] = riwish( prior$nu + N, Phi )
        
        ## update the log likelihood value (for sanity)
        strg$ll[,i] = loglikelihood(elstH, elstCur, N, K, strg$phi0Trc[,i], strg$phi1Trc[,i])
	
        if(i %% 50 == 0){ print( paste("Finished iteration ", i ) ) }
    }
    t2 = Sys.time()
    time = difftime(t2, t1, units="mins")

    ## return output
    return(list(strg=strg, ARu=count_u/nIts, ARr=count_r/nIts, da=da/(nIts*(N-d)), S_u=S_u, S_r=S_r, time=time))
}

plot_mcmc_output <- function(strg, indx, elst){

    K = dim(strg$rTrc)[1] + 1
    N = dim(strg$uTrc)[1]
        
    layout(matrix(c(1,1,2,2,3,1,1,2,2,4,5,6,7,8,9), nrow=3, ncol=5, byrow=TRUE))

    maxv1 = max( strg$uTrc[,1,indx] )
    maxv2 = max( strg$uTrc[,2,indx] )
    minv1 = min( strg$uTrc[,1,indx] )
    minv2 = min( strg$uTrc[,2,indx] )
    cols = rainbow(N)
    cols = sample( cols, N, replace=FALSE)
    plot( NA, xlim=c(minv1, maxv1), ylim=c(minv2, maxv2), pch=20, main="latent coordinates")
    for (i in 1:N){ points( strg$uTrc[i,1,indx], strg$uTrc[i,2,indx], col=cols[i], pch=20 ) }
    umn = apply(strg$uTrc[,,indx], c(1,2), mean)
    text( umn[,1], umn[,2], 1:N)

    plotHypergraph(umn, elst, K)

    plot( colSums(strg$ll)[indx], main="log-likelihood" )

    plot(strg$muTrc[1,indx], strg$muTrc[2,indx], pch=19, main="mu")

    plot(strg$sigTrc[1,1,indx],  type='l', main="sigma[1,1]")
    plot(strg$sigTrc[2,2,indx],  type='l', main="sigma[2,2]")

    plot(strg$rTrc[1,indx], col=1, type='l', ylim=c(0, max(strg$rTrc[,indx])), main="radii")
    for(k in 2:(K-1)){lines(strg$rTrc[k,indx], col=k, type='l')}

    plot(strg$phi0Trc[1,indx], col=1, type='l', ylim=c(0, max(strg$phi0Trc[,indx])), main="phi0")
    for(k in 2:(K-1)){lines(strg$phi0Trc[k,indx], col=k, type='l')}

    plot(strg$phi1Trc[1,indx], col=1, type='l', ylim=c(0, max(strg$phi1Trc[,indx])), main="phi1")
    for(k in 2:(K-1)){lines(strg$phi1Trc[k,indx], col=k, type='l')}
    
}

plot_mcmc_output_us_only <- function(strg, indx, elst){

    K = dim(strg$rTrc)[1] + 1
    N = dim(strg$uTrc)[1]
        
    par(mfrow=c(1,2))

    maxv1 = max( strg$uTrc[,1,indx] )
    maxv2 = max( strg$uTrc[,2,indx] )
    minv1 = min( strg$uTrc[,1,indx] )
    minv2 = min( strg$uTrc[,2,indx] )
    cols = rainbow(N)
    cols = sample( cols, N, replace=FALSE)
    plot( NA, xlim=c(minv1, maxv1), ylim=c(minv2, maxv2), pch=20, main="latent coordinates traceplot", xlab="us[,1]", ylab="us[,2]")
    for (i in 1:N){ points( strg$uTrc[i,1,indx], strg$uTrc[i,2,indx], col=cols[i], pch=20 ) }
    umn = apply(strg$uTrc[,,indx], c(1,2), mean)
    text( umn[,1], umn[,2], 1:N)

    plotHypergraph(umn, elst, K)
    title(main="Posterior mean of latent coordinates")
    
}

plot_mcmc_output_non_us_trcs <- function(strg, indx, elst){

    K = dim(strg$rTrc)[1] + 1
    N = dim(strg$uTrc)[1]

    par(mfrow=c(2,4))
    
    plot( colSums(strg$ll)[indx], type='l', ylab="log-likelihood", xlab="Iteration" )

    plot(strg$phi0Trc[1,indx], col=1, type='l', ylim=c(0, max(strg$phi0Trc[,indx])), ylab="phi0", xlab="Iteration")
    for(k in 2:(K-1)){lines(strg$phi0Trc[k,indx], col=k, type='l')}

    plot(strg$phi1Trc[1,indx], col=1, type='l', ylim=c(0, max(strg$phi1Trc[,indx])), ylab="phi1", xlab="Iteration")
    for(k in 2:(K-1)){lines(strg$phi1Trc[k,indx], col=k, type='l')}

    plot(strg$muTrc[1,indx], strg$muTrc[2,indx], pch=19, xlab="mu[1]", ylab="mu[2]")

    plot(strg$sigTrc[1,1,indx],  type='l', ylab="sigma[1,1]", xlab="Iteration")
    plot(strg$sigTrc[2,2,indx],  type='l', ylab="sigma[2,2]", xlab="Iteration")
    plot(strg$sigTrc[1,2,indx],  type='l', ylab="sigma[1,2]", xlab="Iteration")

    plot(strg$rTrc[1,indx], col=1, type='l', ylim=c(0, max(strg$rTrc[,indx])), ylab="radii", xlab="Iteration")
    for(k in 2:(K-1)){lines(strg$rTrc[k,indx], col=k, type='l')}

    
}

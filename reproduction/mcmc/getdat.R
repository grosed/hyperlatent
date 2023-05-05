## generate hypergraph data 
library("mvtnorm")
source("../mcmc/mcmc_functions.R")
library("graphics")

genHypergraph <- function( N, K, r, phi0, phi1, mu, sigma, add_noise = TRUE, add_noise_full = FALSE, bkstn=FALSE )
{
    d = length(mu) # latent dimension
    if( d > 3 ){stop("Can only gen graph for d = 2 or 3")}
    if( d==3 & bkstn==TRUE ){stop("Bookstein not coded for d=3, only for d=2")}

    ## generate coordinates (only for d=2)
    us = rmvnorm(N, mean=mu, sigma=sigma)
    if (bkstn == TRUE)
    {
        ## get bookstein coordinates
        bk = getBookstein(us)
        us = bk$u

        ## update mu and sigma
        mu = bk$c * bk$R %*% (mu - bk$b)
        sigma = (bk$c^2) * bk$R %*% sigma %*% t(bk$R)
    }

    ## induce graph
    edges = list() # initialise list for storing edges
    ## loop over edge sizes and build up graph
    for (k in 2:K) 
    {
        if ( d == 2){
            cech_tmp = calc_cech_R2(rmax=1.1*r[k-1], t(us), N, 2, k-1)
        } else if ( d == 3 ){
            cech_tmp = calc_cech_R3(rmax=1.1*r[k-1], t(us), N, 3, k-1)
        }
        Ilenk = lapply(cech_tmp$nodes, length) == k
        Ifilt = unlist(cech_tmp$filts) <= r[k-1]
        if (sum( Ilenk & Ifilt ) > 0 )
        {
            emat = matrix(unlist(cech_tmp$nodes[ Ilenk & Ifilt]), ncol=k, byrow=TRUE) + 1
            elst = lapply( seq_len(nrow(emat)), function(i) emat[i,] )

            ## update edges
            edges = append( edges, elst )
        }
        
    }

    edges = lapply( edges, sort) # make sure i<j<... etc
    ## add noise
    if (add_noise == TRUE )
    {
        if (add_noise_full==TRUE){
            edges = addnoise_full( edges, N, K, phi0, phi1 )
        } else {
            edges = addnoise( edges, N, K, phi0, phi1 )
        }
    }
    edges = lapply( edges, sort) # make sure i<j<... etc
    edges = unique(edges)
    
    return( list(us=us, mu=mu, sigma=sigma, elst=edges) )
}

genSimplicialHypergraph <- function( N, K, r, phi0, phi1, mu, sigma, add_noise = TRUE, add_noise_full = FALSE, bkstn=FALSE )
{
    d = length(mu) # latent dimension
    if( d > 3 ){stop("Can only gen graph for d = 2 or 3")}
    if( d==3 & bkstn==TRUE ){stop("Bookstein not coded for d=3, only for d=2")}

    ## generate coordinates (only for d=2)
    us = rmvnorm(N, mean=mu, sigma=sigma)
    if (bkstn == TRUE)
    {
        ## get bookstein coordinates
        bk = getBookstein(us)
        us = bk$u

        ## update mu and sigma
        mu = bk$c * bk$R %*% (mu - bk$b)
        sigma = (bk$c^2) * bk$R %*% sigma %*% t(bk$R)
    }

    ## induce graph
    if ( d == 2){
        cech_tmp = calc_cech_R2(rmax=1.1*r, t(us), N, 2, K-1)
    } else if ( d == 3 ){
        cech_tmp = calc_cech_R3(rmax=1.1*r, t(us), N, 3, K-1)
    }
    Ifilt = unlist(cech_tmp$filts) <= r
    Ilen = unlist(lapply(cech_tmp$nodes,length)) > 1
    if (sum( Ilen & Ifilt ) > 0 ){edges = cech_tmp$nodes[ Ifilt & Ilen ]
    } else { edges=list() }
    edges = lapply( edges, function(x){ x + 1 } )
    edges = lapply( edges, sort) # make sure i<j<... etc
    ## add noise
    if (add_noise == TRUE )
    {
        if (add_noise_full==TRUE){
            edges = addnoise_full( edges, N, K, phi0, phi1 )
        } else {
            edges = addnoise( edges, N, K, phi0, phi1 )
        }
    }
    edges = lapply( edges, sort) # make sure i<j<... etc
    edges = unique(edges)
    
    return( list(us=us, mu=mu, sigma=sigma, elst=edges) )
}

## plot hypergraph functions
plotHypergraph <- function(us, elst, K)
{

    N = dim(us)[1]
    grphks = unlist( lapply( elst, length ) ) # edge dimensions
    
    ## plot points
    rng = c( min(us), max(us) )
    plot( us[,1], us[,2], pch=16, cex=1.5, xlim=rng, ylim=rng)

    ## plot ijk
    e_ijk = elst[ grphks==3 ]
    n_ijk = length(e_ijk)

    if (n_ijk !=0 )
    {
        for ( i in 1:n_ijk )
        {
            ecur = e_ijk[[i]]
            polygon( us[ecur, 1], us[ecur, 2], col="skyblue", border="skyblue", lwd=2, density=30, fillOddEven=TRUE )
        }
    }
    
    ## plot ij
    e_ij = elst[ grphks==2 ]
    n_ij = length(e_ij)

    if (n_ij !=0 )
    {
        for( i in 1:n_ij )
        {
            ecur = e_ij[[i]]
            segments(us[ ecur[1] ,1], us[ ecur[1] ,2], us[ ecur[2] ,1], us[ ecur[2] ,2], lwd=2, col="darkorange")
        }
    }

    ## add labels
    points( us[,1], us[,2], pch=16, cex=1.5 ) #for ordering
    text( us[,1], us[,2], labels=seq(1,N), cex=1.25, pos=3 )
    
}

## simulate new connections for a hypergraph
genHyperG_with_NewConnection <- function(Nstr, prms, periph=FALSE, addnoise=TRUE)
{
                                        # Nstr is number of new nodes
                                        # prms is list of model parameters (ie posterior samples)

    ## initial points have index 1:prms$N
    ## new points have index (prms$N+1):(prms$N+Nstr)
    d = length( prms$mu )
    
    ## sample new points
    if (periph==FALSE){
        uStr =  rmvnorm(Nstr, mean=prms$mu, sigma=prms$sigma)
        uFull = rbind( prms$us, uStr ) # use these points to simulate new graph connections
    } else if (periph==TRUE){
        Ntot = prms$N + Nstr # sample all nodes and take Nstr furthest pts
        uStr = rmvnorm(Ntot, mean=prms$mu, sigma=prms$sigma)
        ## calculate and sort distances from mu
        uptssrt = sort(apply( uStr, 1, function(x){ dist( rbind(x, prms$mu) ) } ), index.return=TRUE)
        ## take 'new' us to be those furthest from mean mu
        uStr = uStr[uptssrt$ix[(prms$N+1):Ntot], ]
        uFull = rbind( prms$us, uStr )
    }
    NFull = prms$N + Nstr
    
    ## induce graph
    edges = list() # initialise list for storing edges
    ## loop over edge sizes and build up graph
    for (k in 2:prms$K) 
    {
        if ( d == 2){
            cech_tmp = calc_cech_R2(rmax=1.1*prms$r[k-1], t(uFull), NFull, 2, k-1)
        } else if ( d == 3 ){
            cech_tmp = calc_cech_R3(rmax=1.1*prms$r[k-1], t(uFull), NFull, 3, k-1)
        }
        Ilenk = lapply(cech_tmp$nodes, length) == k
        Ifilt = unlist(cech_tmp$filts) <= prms$r[k-1]
        if (sum( Ilenk & Ifilt ) > 0 )
        {
            emat = matrix(unlist(cech_tmp$nodes[ Ilenk & Ifilt]), ncol=k, byrow=TRUE) + 1
            elst = lapply( seq_len(nrow(emat)), function(i) emat[i,] )
            ## update edges
            edges = append( edges, elst )
        }
    }

    ## add noise
    if (addnoise==TRUE){ edges = addnoise( edges, NFull, prms$K, prms$phi0, prms$phi1 ) }

    ## sanity check hyperedges are unique
    edges = lapply(edges, sort)
    edges = unique(edges)
    
    ## subset the hyperedges from full hypergraph containing the new indexes
    idNew = (prms$N+1):NFull
    nedge = length(edges)

    if (nedge > 0){ # if hypergraph not empty

        ## get hypergraph containing 'old' nodes only
        idOld = 1:prms$N
        idOldNode = unlist(lapply( edges, function(i){ all(i %in% idOld) } ))
        edges_old = edges[idOldNode]
        
        ## subset hyperedges contain new index 
        idNewNode = unlist(lapply( edges, function(i){ any(idNew %in% i) } ))
        
        ## separate hyperedges into only new, and between new and old
        if ( sum(idNewNode) > 0 ){
            etmp = edges[idNewNode]
            
            nNewNodes = lapply( etmp, function(i){ sum(idNew %in% i) } )
            lEdge = unlist( lapply( etmp, length ) ) # length of edges
            edges_new = etmp[ lEdge == nNewNodes ]

            ## 2) between new and old
            edges_btwn = etmp[ lEdge != nNewNodes ]
        } else {
            edges_new = list()
            edges_btwn = list()
        }

    } else { # if hypergraph empty
        edges = list()
        edges_new = list()
        edges_btwn = list()
    }
    
    return( list(eOld = edges_old, eNew = edges_new, eBtwn = edges_btwn, us=uFull) )
}

genHyperG_given_prms <- function(prms, addnoise=TRUE)
{
                                        # prms is list of model parameters (ie posterior samples)

    ## initial points have index 1:prms$N
    ## new points have index (prms$N+1):(prms$N+Nstr)
    d = length( prms$mu )
       
    ## induce graph
    edges = list() # initialise list for storing edges
    ## loop over edge sizes and build up graph
    for (k in 2:prms$K) 
    {
        if ( d == 2){
            cech_tmp = calc_cech_R2(rmax=1.1*prms$r[k-1], t(prms$us), prms$N, 2, k-1)
        } else if ( d == 3 ){
            cech_tmp = calc_cech_R3(rmax=1.1*prms$r[k-1], t(prms$us), prms$N, 3, k-1)
        }
        Ilenk = lapply(cech_tmp$nodes, length) == k
        Ifilt = unlist(cech_tmp$filts) <= prms$r[k-1]
        if (sum( Ilenk & Ifilt ) > 0 )
        {
            emat = matrix(unlist(cech_tmp$nodes[ Ilenk & Ifilt]), ncol=k, byrow=TRUE) + 1
            elst = lapply( seq_len(nrow(emat)), function(i) emat[i,] )
            ## update edges
            edges = append( edges, elst )
        }
    }

    ## add noise
    if (addnoise==TRUE){ edges = addnoise( edges, prms$N, prms$K, prms$phi0, prms$phi1 ) }

    ## sanity check hyperedges are unique
    edges = lapply(edges, sort)
    edges = unique(edges)
    
    return( list(elst=edges, us=prms$us) )
}

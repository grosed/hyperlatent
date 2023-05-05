## code to simulate hypergraphs for each type of misspec
library("mvtnorm")
source("../mcmc/mcmc_functions.R") 
library("igraph")

## 1) function for u misspecification
## 2) function for non-homogeneous errors 
## 3) function for non-homogeneous radii
## 4) function for simplicial hypergraph

######################### latent misspec function ####################

getlatentcoords <- function(N, d, prms, ucase)
{
                                        # determine latent coordinates
                                        # N is number of nodes
                                        # d is latent dimension
                                        # prms contains distribution prms
                                        # ucase is 1) "nomiss" - from Normal
                                        # 2) "clust" - mixture of normalss
                                        # 3) "unif" - uniform distn

    if ( ucase %in% c("nomiss", "clust", "unif") == FALSE )
    {
        stop("ucase is not valid. Can be 'nomiss', 'clust' or 'unif'")
    }

    if (ucase == "nomiss"){

        ## sample coordinates from N(mu, sigma)
        us = rmvnorm(N, mean=prms$mu, sigma=prms$sigma)
        
    } else if (ucase == "clust") {

        ## sample coordinates from mixture of normals
        us = matrix( NA, nrow=N, ncol=d )
        C = length(prms$lvec) # number of clusters
        id = sample(1:C, prob=prms$lvec, size=N, replace=TRUE)

        for (n in 1:N) # loop ok as one off
        {
            us[n,] = rmvnorm(1, mean=prms$mu[,id[n]], sigma=prms$sigma[,,id[n]] )
        }

        
    } else if (ucase == "unif" ) {

        ## sample coordinates uniformly on [0,1]^d
        us = matrix(runif(N*d), nrow=N, ncol=d)
        
    }

    return( us )
}

####################### non-homog error function ####################

addnonhomnoise <- function( edges, N, K, prms, phicase )
{

    elen = unlist( lapply( edges, length) )
    eout = list()
    
    ## phicase is 'nomiss', 'nonhom'
    if ( phicase == "nomiss" )
    {

        for (k in 2:K)
        {
            ## separate edges into those in and not in the graph
            ek = edges[elen==k]
            efullk = combn( N, k )
            efullk = lapply(seq_len(ncol(efullk)), function(i) efullk[,i]) # list of all e_k
            efullk = lapply( efullk, sort)
            ek = lapply( ek, sort)

            ids_ek = lapply( ek, function(x){which(unlist(lapply(efullk, identical, as.integer(x))))} )
            abs_ek = efullk[ - unlist(ids_ek) ]            
                        
            ## perturb edges in g, ek
            nek = length(ek)
            idek = rbinom( nek, 1, prms$phi1[k-1] )

            ## perturb edges not in g, abs_ek
            nabs_ek = length(abs_ek)
            idabsek = rbinom( nabs_ek, 1, prms$phi0[k-1] )

            ## take all edges present in the graph
            ek = append(ek[idek==0], abs_ek[idabsek==1])
            eout = append( eout, ek )
        }       
        
    } else if ( phicase == "nonhom" )
    {
        
        for (k in 2:K)
        {
            ## separate edges into those in and not in the graph
            ek = edges[elen==k]
            efullk = combn( N, k )
            efullk = lapply(seq_len(ncol(efullk)), function(i) efullk[,i]) # list of all e_k
            
            ids_ek = lapply( ek, function(x){which(unlist(lapply(efullk, identical, as.integer(x))))} )
            abs_ek = efullk[ - unlist(ids_ek) ]
            
            ## perturb edges in g, ek
            nek = length(ek)
            idek = rep(NA, nek)
            for (l in 1:nek)
            {
                ptmp = rbeta(1, prms$a1vec[k-1], prms$b1vec[k-1])
                idek[l] = rbinom(1, 1, ptmp)
            }

            ## perturb edges not in g, abs_ek
            nabs_ek = length(abs_ek)
            idabsek = rep(NA, nabs_ek)
            for (l in 1:nabs_ek)
            {
                ptmp = rbeta(1, prms$a0vec[k-1], prms$b0vec[k-1])
                idabsek[l] = rbinom(1, 1, ptmp)
            }

            ## take all edges present in the graph
            ek = append(ek[idek==0], abs_ek[idabsek==1])
            eout = append( eout, ek )
        }
        
    } else { stop("Invalid choice for phicase. Can be 'nomiss' or 'nonhom' ") }

    return( eout )
}

######################## gen hypergraph shell #######################

getk2edge <- function(N, us, rmat, k, edges)
{
                                        # N is number of nodes
                                        # rmat contains radii (dim Nx(K-1)
                                        # k is current e order
                                        # edges is current edge list
    k2comb = combn( N, 2 )
    ncomb = dim(k2comb)[2]
    for (i in 1:ncomb)
    {
        dtmp = dist( us[k2comb[,i],] )
        if (dtmp <= rmat[k2comb[1,i],k-1] + rmat[k2comb[2,i],k-1] )
        {
            edges = append( edges, list(k2comb[,i]) )
        }
    }

    return(edges)
}

## calculates points of intersection for two disks different radii
getintpts <- function(us, rmat, k, id)
{
                                        # find points of intersection
                                        # for circles specified by id

    if (length(id) != 2){stop("function is for 2 circles only")}

    ##  get coordinates and radii
    xa = us[id[1],1]
    ya = us[id[1],2]
    xb = us[id[2],1]
    yb = us[id[2],2]
    ra = rmat[id[1],k-1]
    rb = rmat[id[2],k-1]

    ## calculate y intersection point
    const = ( (ra^2 -rb^2) - (xa^2-xb^2) - (ya^2-yb^2) ) / (2*(xb-xa))
    A = ((ya-yb)/(xb-xa))^2 + 1
    B = -2*(((ya-yb)/(xb-xa))*(xa-const) + ya)
    C = ya^2 + (xa - const)^2 - ra^2

    if ( (B^2 - 4*A*C) != 0 ){ # ie. two solns for y
        
        yst1 = (-B + sqrt(B^2 - (4*A*C)))/(2*A)
        yst2 = (-B - sqrt(B^2 - (4*A*C)))/(2*A)
        
        ## calculate x intersection point
        xst1 = const + yst1 * ((ya-yb)/(xb-xa))
        xst2 = const + yst2 * ((ya-yb)/(xb-xa))
        
    } else { # 1 soln for y, solve for x first

        ## calculate x intersection point
        const = ( (ra^2 -rb^2) - (xa^2-xb^2) - (ya^2-yb^2) ) / (2*(yb-ya))
        A = ((xa-xb)/(yb-ya))^2 + 1
        B = -2*(((xa-xb)/(yb-ya))*(ya-const) + xa)
        C = xa^2 + (ya - const)^2 - ra^2

        xst1 = (-B + sqrt(B^2 - (4*A*C)))/(2*A)
        xst2 = (-B - sqrt(B^2 - (4*A*C)))/(2*A)
        
        ## calculate y intersection point
        yst1 = const + xst1 * ((ya-yb)/(xb-xa))
        yst2 = const + xst2 * ((ya-yb)/(xb-xa))

    }
        
    ## return answers
    return(list(int1=c(xst1, yst1), int2=c(xst2, yst2)))
}

genHyp_shell <- function( N, K, d, prms, ucase, rcase, phicase )
{
                                        # N is number of nodes
                                        # K is max hyperedge order
                                        # prms contains model parameters
                                        # d is latent dimension
                                        # ucase can be "nomiss", "clust", "unif"
                                        # rcase can be "nomiss", "nonhom", "simp"
                                        # phicase can be "nomiss", "nonhom"
    
    if( d > 3 ){stop("Can only gen graph for d = 2 or 3")}

    ## generate coordinates (only for d=2)
    us = getlatentcoords(N, d, prms, ucase)

    ## determine the hypergraph from U and r
    if ( rcase %in% c( "nomiss", "simp" ) == TRUE )
    {

        ## find edges as in genhyp function
        ## in simplicial case, r is a constant vector of length K-1
        
        ## induce graph
        edges = list() # initialise list for storing edges
        ## loop over edge sizes and build up graph
        for (k in 2:K) 
        {
            if ( d == 2){
                cech_tmp = calc_cech_R2(rmax=prms$r[k-1], t(us), N, 2, K-1)
            } else if ( d == 3 ){
                cech_tmp = calc_cech_R3(rmax=prms$r[k-1], t(us), N, 3, K-1)
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

        edges = lapply( edges, sort)
        
    } else if ( rcase == "nonhom" ) { 

        ## populate rmat from priors
        rmat = matrix(NA, ncol=K-1, nrow=N)
        for (k in 2:K)
        {
            rmat[,k-1] = rbeta(N, prms$ravec[k-1], prms$rbvec[k-1]) * (prms$rlim[k]-prms$rlim[k-1]) + prms$rlim[k-1] 
        }

        print( max(rmat[,1] ) )
        print( max(rmat[,2] ) )
              
              
        if( K > 3 ){stop("No code for k order intersection of disks with different radii when k>3")}
        
        ## in this setting the radii are different for each nodes
        ## i=1,2,...N nodes, and for each k=2,3,...,K radii

        ## determine pairwise edges
        edges = list() # empty list for edges
        edges = getk2edge(N, us, rmat, 2, edges) # k=2 or 3 for when radii unique to each node

        ## determine order 3 edges (if included in graph)
        if ( K == 3 )
        {
            etmp = list()
            etmp = getk2edge(N, us, rmat, 3, etmp) 
            etmp = matrix(unlist( etmp ), ncol=2, byrow=TRUE )
            ## take all triples of these edges
            gtmp = graph_from_edgelist( etmp, directed=FALSE )
            tris = t(apply(matrix( triangles(gtmp), nrow=3 ),2,sort))
            ntris = dim(tris)[1]

            ## now check which order 3 edges are present
            for (tr in 1:ntris)
            {

                ## only find intersection in 1st case
                ## each pair intersects
                
                ## check no disks are contained within eachother
                d12 = dist( us[tris[tr, c(1,2)],] )
                d13 = dist( us[tris[tr, c(1,3)],] )
                d23 = dist( us[tris[tr, c(2,3)],] )
                
                if ( d12 < (max(rmat[tris[tr,c(1,2)],2]) - min(rmat[tris[tr,c(1,2)],2])) | d13 < (max(rmat[tris[tr,c(1,3)],2]) - min(rmat[tris[tr,c(1,3)],2])) | d23 < (max(rmat[tris[tr,c(2,3)],2]) - min(rmat[tris[tr,c(1,2)],2])) )
                {
                    
                    edges = append(edges, list(tris[tr,]))
                    
                } else {
                    
                    ## find intersection pts of 1st two circles
                    inttmp = getintpts(us, rmat, 3, tris[tr,c(1,2)])

                    ## check distance between int and 3rd disk center
                    d1 = dist( matrix( c(inttmp$int1, us[tris[tr,3],]), ncol=2, byrow=TRUE ) )
                    d2 = dist( matrix( c(inttmp$int2, us[tris[tr,3],]), ncol=2, byrow=TRUE ) )

                    ## add hyperedge if dist < 3rd disk radius
                    if ( d1 <= rmat[tris[tr,3],2] | d2 <= rmat[tris[tr,3],2] )
                    {
                        edges = append(edges, list(tris[tr,]))
                    }
                    
                }
             
            }

        }
        edges = lapply( edges, sort)
    } else { stop("Invalid choices for rcase. Can be 'nomiss', 'nonhom' or 'simp'") }
    
    ## add noise to hyperedges
    edges = addnonhomnoise( edges, N, K, prms, phicase )
    
    return( list(us=us, elst=edges) )
}

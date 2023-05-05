#library(Rcpp)
#Sys.setenv("PKG_CXXFLAGS"="-std=c++11") 
#Sys.setenv("PKG_CXXFLAGS"="-I/usr/shared_apps/packages/anaconda3-2020.11/include")
#sourceCpp('cech_complex.cpp') 
#sourceCpp('miniball.cpp')
library('igraph')
library('geometry') # for 3d convex hull, bk calculation
library('dplyr')

convert_l_to_df <- function(lst, k){
	if (length(lst) > 0){
	   df = do.call(rbind.data.frame, lst)
	   names(df) = NA
	} else {
	   df = data.frame( matrix( vector(), 0, k ) )
	   names(df) = NA
	}
	return(df)
}		

###################### INDUCE GRAPH ########################################

get_graph_from_cech <- function(us, r)
{
                                        # calculate a graph using the cech complex
                                        # us as N x 2 matrix of latent coordinates
                                        # r is vector of radii

    K = length(r) + 1 # (don't have r_1)
    N = dim(us)[1]
    latd = dim(us)[2]

    edges = list() # initialise list for storing edges
    ## loop over edge sizes and build up graph
    for (k in 2:K) 
    {
        ##print( r[k-1] )
        if( latd==2 ){
            cech_tmp = calc_cech_R2(rmax=r[k-1], t(us), N, 2, k-1)
        } else if (latd==3){
            cech_tmp = calc_cech_R3(rmax=r[k-1], t(us), N, 3, k-1)
        }
        Ilenk = lengths(cech_tmp$nodes) == k
        Ifilt = unlist(cech_tmp$filts) <= r[k-1]
        if (sum( Ilenk & Ifilt ) > 0 )
            {
                emat = matrix(unlist(cech_tmp$nodes[ Ilenk & Ifilt]), ncol=k, byrow=TRUE) + 1
                elst = lapply( seq_len(nrow(emat)), function(i) emat[i,] )

                ## update edges
                edges = append( edges, elst )
            }        
    }
    edges = lapply( edges, sort ) # sort hyperedges into ascending order
    
    return( edges )
}

get_simplicial_graph_from_cech <- function(us, r, K)
{
                                        # calculate a graph using the cech complex
                                        # us as N x 2 matrix of latent coordinates
                                        # r is vector of radii

    N = dim(us)[1]
    latd = dim(us)[2]

    ## since simplicial, only one cech evaluation
    if( latd==2 ){
        cech_tmp = calc_cech_R2(rmax=1.1*r, t(us), N, 2, K-1)
    } else if (latd==3){
        cech_tmp = calc_cech_R3(rmax=1.1*r, t(us), N, 3, K-1)
    }
    Ilen = lapply(cech_tmp$nodes, length) > 1 # remove hyperdges with one node
    Ifilt = unlist(cech_tmp$filts) <= r
    if (sum( Ilen & Ifilt ) > 0 ){
        edges = cech_tmp$nodes[ Ilen & Ifilt ]
        edges = lapply(edges, function(x){x+1})
    } else {
        edges = list()
    }    
    edges = lapply( edges, sort ) # sort hyperedges into ascending order
    
    return( edges )
}

get_graph_from_cech_order_k <- function(us, r, k)
{
                                        # calculate a graph using the cech complex
                                        # us as N x 2 matrix of latent coordinates
                                        # r is vector of radii
                                        # k is 'current' hyperedge order in (2,3,...,K)
    N = dim(us)[1]
    latd = dim(us)[2]

    if( latd==2 ){
        cech_tmp = calc_cech_R2(rmax=r[k-1], t(us), N, 2, k-1) 
    } else if (latd==3){
        cech_tmp = calc_cech_R3(rmax=r[k-1], t(us), N, 3, k-1)
    }
    Ilenk = lengths(cech_tmp$nodes) == k
    Ifilt = unlist(cech_tmp$filts) <= r[k-1]
    if (sum( Ilenk & Ifilt ) > 0 )
    {
        emat = matrix(unlist(cech_tmp$nodes[ Ilenk & Ifilt]), ncol=k, byrow=TRUE) + 1
        edges = lapply( seq_len(nrow(emat)), function(i) emat[i,] )
    } else {
        edges = list()
    }

    edges = lapply( edges, sort ) # sort hyperedges into ascending order
        
    return( edges )
}

sample_ek_ind <- function(N, k){
    ## sample random hyperedge of order k
    ## since unique, can sample k from 1:N without replacement and sort
    ind = sort(sample(1:N, k, replace=FALSE))
}

pet_edges <- function( elst, N, nflp, k )
{
                                        # elst - list of k order edges
                                        # N is number of nodes
                                        # nflp is number of e to flip
                                       
    elstadd = list()
    
    ## flip edges
    if (nflp > 0)
    {
        for (i in 1:nflp) # loop over edges
        {
            ## sample an index (sample until we have a new index)
            repeat
            {
                I = sample_ek_ind( N, k )
                if (list(I) %in% elst == FALSE){ break }
            }
            
            ## add in edge
            elstadd = append( elstadd, list(I) )
        }
    }

    return( elstadd ) # return perturbed matrix
}

addnoise <- function(edgeG, N, K, phi0, phi1)
{
    ## function to perturb edges in the graph
    edgeG = lapply(edgeG, sort) # just to be sure
    
    elstOut = list() # new list of edges
    elen = unlist( lapply( edgeG, length ) ) # get edge lengths
    
    for (k in 2:K)
    {
        ek = edgeG[ elen == k] # extract edges of current size
        nk = length(ek) # number of order k edges
        
        ## add phi0 noise (ie include additional hyperedges)
        nflp = choose(N,k) - nk #number I can flip from (remove 1's)
        nflp = rbinom(1, nflp, phi0[k-1]) # number to flip
        e_add = pet_edges( ek, N, nflp, k ) # e's to add
        elstOut = append(elstOut, e_add) # add new edges

        ## add phi1 noise (ie remove existing hyperedges)
        flp = rbinom(nk, 1, phi1[k-1])
        elstOut = append(elstOut, ek[ flp==0 ] ) # if 0, do not change
    }

    elstOut = unique( elstOut )
    elstOut = lapply( elstOut, sort) 

    return( elstOut )
}

addnoise_full<- function(edgeG, N, K, phi0, phi1)
{
    ## function to perturb edges in the graph
    edgeG = lapply(edgeG, sort) # just to be sure    
    elstOut = list() # new list of edges
    elen = unlist( lapply( edgeG, length ) ) # get edge lengths
    for (k in 2:K)
    {
        ek = edgeG[ elen == k] # extract edges of current size
        nk = length(ek) # number of order k edges
        
        if (nk > 0){
            ## all edges
            efull = combn(N, k)
            ntot = dim(efull)[2]
            idin = rep(NA, ntot)
            ## ek not in graph
            for (i in 1:ntot){
                ith_id = lapply(1:nk, function(j){ sum(efull[,i] %in% ek[[j]] )})
                idin[i]= sum( unlist(ith_id) == k )
            }
            e_abs = matrix(efull[,idin==0], ncol = sum(idin==0) ) # edges not in g
            e_in = matrix(efull[,idin==1], ncol = sum(idin==1) ) # edges in g
            
            ## add phi0 noise (include additional hyperedges)
            idflp = rbinom( dim(e_abs)[2], 1, phi0[k-1] )
            e_add = e_abs[, idflp==1]
            if( is.numeric( ncol(e_add)) == TRUE ){
            e_add = lapply(seq_len(ncol(e_add)), function(i) e_add[,i]) } else {e_add = list() }
            elstOut = append(elstOut, e_add) # add new edges

            ## add phi1 noise (remove present hyperedges)
            idflp = rbinom( dim(e_in)[2], 1, phi1[k-1] )
            e_inc = e_in[, idflp==0]
            if( is.numeric( ncol(e_inc)) == TRUE ){
            e_inc = lapply(seq_len(ncol(e_inc)), function(i) e_inc[,i]) } else { e_inc = list() }

            elstOut = append(elstOut, e_inc) # if idflp=0, keep hyperedge
        } else {
            ## no edges in graph, can only add edges into g
            e_abs =  combn(N, k)
            nabs = dim(e_abs)[2]
            idflp = rbinom( nabs, 1, phi0[k-1])
            e_add = e_abs[, idflp==1 ]
            if( is.numeric(ncol(e_add)) == TRUE ){
            e_add = lapply(seq_len(ncol(e_add)),function(i) e_add[,i])
            } else {e_add = list() }
            
            elstOut = e_add
        }
    }

    elstOut = lapply( elstOut, sort) # just to be sure

    return( elstOut )
}

###################### LOG LIKELIHOOD #####################################x###

loglikelihood <- function(edgeH, edgeG, N, K,  phi0, phi1)
{
                                        # calculate log likelihood from edge lists
                                        # edgeH is observed graph
                                        # edgeG is current graph
                                        # N is number of nodes in graph
                                        # K is max edge order
                                        # phi0 is model noise (K-1)
                                        # phi1 is model noise (K-1)

    ll = rep(0, K-1)

    ## create list of edge lengths
    lenH = lengths( edgeH )
    lenG = lengths( edgeG )
    
    ## loop over edge orders
    for (k in 2:K)
    {
        nCmb = choose(N, k) # number of possible edges of this order

        ## need to extract edges of specific order from each list
        Hktmp = edgeH[ lenH == k ]
        Gktmp = edgeG[ lenG == k ]

        ## find edge counts
        M11 = dim(intersect( convert_l_to_df(Gktmp, k), convert_l_to_df(Hktmp, k)))[1]
        M10 = length( Gktmp ) - M11 # edges in induced graph, not in observed graph
        M01 = length( Hktmp ) - M11 # edge in observed graph, not in induced graph
        M00 = nCmb - (M11 + M10 + M01) # edges is neither graph

        ## update log likelihood
        ll[k-1] = (M10 * log(phi1[k-1])) + (M11 * log(1-phi1[k-1])) + (M01 * log(phi0[k-1])) + (M00 * log(1-phi0[k-1] ))
    }
    
    return( ll )
}

ll_order_k <- function(edgeH, edgeG, N, k,  phi0, phi1)
{
                                        # calculate log likelihood from edge lists
                                        # edgeH is observed graph
                                        # edgeG is current graph
                                        # N is number of nodes in graph
                                        # k is 'current' hyperedge order in (2,3,...,K)
                                        # phi0 is model noise (K-1)
                                        # phi1 is model noise (K-1)
    
    ## create list of edge lengths
    lenH = lengths( edgeH )
    lenG = lengths( edgeG )
    
    nCmb = choose(N, k) # number of possible edges of this order

    ## extract edges of order k from each list
    Hktmp = edgeH[ lenH == k ]
    Gktmp = edgeG[ lenG == k ]

    ## find edge counts
    M11 = dim(intersect( convert_l_to_df(Gktmp, k), convert_l_to_df(Hktmp, k)))[1]
    M10 = length( Gktmp ) - M11 # edges in induced graph, not in observed graph
    M01 = length( Hktmp ) - M11 # edge in observed graph, not in induced graph
    M00 = nCmb - (M11 + M10 + M01) # edges is neither graph

    ## update log likelihood
    ll = (M10 * log(phi1[k-1])) + (M11 * log(1-phi1[k-1])) + (M01 * log(phi0[k-1])) + (M00 * log(1-phi0[k-1] ))
    
    return( ll )
}


###################### BOOKSTEIN COORDINATES ########################################

eucdist <- function(a,b) {
    d = sqrt( sum( (a-b)^2 ) )
}

getBookstein <- function(u, d=2, PCA=FALSE, ids=NA)
{
                                        # function to calculate bookstein coordinates
                                        # u is Nx2 array of coordinates
                                        # d is the dimension of the latent space, can be 2 or 3
                                        # PCA = FALSE -> fixed points found on convex hull
                                        # PCA = TRUE -> fixed points found on principal axis
                                        # if ids == NA, choose given pca or convex hull

    if (d==2){
        N = nrow(u)
        ## select points to be anchored, either by convex hull of pt cloud or PCA
        if (any(is.na(ids))==FALSE){
            IB = ids
        } else {
            if ( PCA == FALSE ){
                ch = chull(u) # convex hull of point cloud
                IB = ch[1:2] # choose anchors from this point cloud
            } else if ( PCA == TRUE ){
                pca = prcomp(u) # apply PCA to coordinates
                isml = which.min(pca$sdev) # choose smallest axis of variation
                evec = pca$rotation[,isml] # corrosponding eigen vector
                pt1 = 2*pca$sdev[isml]*evec + pca$center # end point of axis of variation line
                pt2 = (-2)*pca$sdev[isml]*evec + pca$center # other end point
                ## calculate distances from points to line
                ds1 = apply(u, 1, function(i) eucdist(a=i, b=pt1))
                ds2 = apply(u, 1, function(i) eucdist(a=i, b=pt2))
                IB1 = (sort( ds1, decreasing=FALSE, index.return=TRUE)$ix)[1] # v closest to to pt1
                IB2 = (sort( ds2, decreasing=FALSE, index.return=TRUE)$ix)[1] # v closest to to pt2
		if (IB1 == IB2){ IB2 = (sort( ds2, decreasing=FALSE, index.return=TRUE)$ix)[2] }
                IB = c(IB1, IB2)
            }
        }
        IBsrt = (sort(u[IB,1], index.return=TRUE))$ix
        IB = IB[IBsrt] # make sure points are ordered according to x (prevents reflection in y axis)
        
        ## calculate bookstein
        uBook = matrix( NA, N, 2)

        ## define anchors
        uBook[IB[1],] = c(-.5, 0) 
        uBook[IB[2],] = c(.5, 0)

        ## transform other points
        xdiff = u[IB[2],1] - u[IB[1],1]
        ydiff = u[IB[2],2] - u[IB[1],2]
        D12 = (xdiff)^2 + (ydiff)^2
        const = 1/sqrt(D12)
        angle = atan( ydiff / xdiff )
        R = matrix( c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol=2, byrow=FALSE )
        b = 0.5*c( sum(u[IB[1:2],1]), sum(u[IB[1:2],2]) ) 

        for (i in (1:N)[-IB])
        {
            uBook[i,] = const * R %*% ( u[i,] - b )
        }

        ## return output (including rescaling matrices)
        if (any(is.na(ids))==TRUE){
            if (PCA == FALSE){
                return( list( u=uBook, id=IB, c=const, R=R, b=b) )
            } else if (PCA == TRUE){
                return( list( u=uBook, id=IB, pca=pca, c=const, R=R, b=b) ) 
            }
        } else {
            return( list( u=uBook, id=IB, c=const, R=R, b=b) )
        }

    } else if (d==3){

        N = nrow(u)
        ## sample the anchor points according to point cloud
        ## function works for an N x d matrix
        ch = convhulln(u) # convex hull of point cloud
        IB = ch[1,] # index of anchors

        ## apply bk transformation
        xIB1 = u[IB[1],]#c(-0.5,0,0)
        xIB2 = u[IB[2],]#c(0.5,0,0)
        ## step 1: calculate ws
        ws = matrix( NA, nrow=N, ncol=3 )
        for (i in 1:N){
            ws[i,] = u[i,] - (xIB1 + xIB2)/2
        }

        D12 = 2 * sqrt( sum(ws[IB[2],]^2) )

        ## calculate angles
        theta = atan( ws[IB[2],2] / ws[IB[2],1] )
        omega = atan( ws[IB[2],3] / sqrt( ws[IB[2],1]^2 + ws[IB[2],2]^2 ) )
        phinum = (ws[IB[2],1]^2 + ws[IB[2],2]^2)*ws[IB[3],3] - ( (ws[IB[2],1]*ws[IB[3],1]) + (ws[IB[2],2]*ws[IB[3],2]) )*ws[IB[2],3]
        phidenom = sqrt( sum(ws[IB[2],]^2) ) * ( (ws[IB[2],1]*ws[IB[3],2]) - (ws[IB[3],1]*ws[IB[2],2]) )
        phi = atan( phinum / phidenom )

        ## calculate rotation matrices
        Rx = matrix(c(1,0,0,0,cos(phi),-sin(phi),0,sin(phi),cos(phi)), ncol=3)
        Ry = matrix(c(cos(omega),0,-sin(omega),0,1,0,sin(omega),0,cos(omega)), ncol=3)
        Rz = matrix(c(cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1), ncol=3)

        ## find bookstein:
        uBook = matrix(NA, nrow=N, ncol=3)
        uBook[(1:N)[-IB[1:2]],] = t(Rx %*% Ry %*% Rz %*% t(ws[(1:N)[-IB[1:2]],])) / D12
        uBook[IB[1],] = c(-0.5,0,0)
        uBook[IB[2],] = c(0.5,0,0)
        
        return( list( u=uBook, id=IB ) )
        
    } else {
        stop("Bookstein coordinates only coded for d=2 and d=3")
    }

}

############################## GMDS init for coordinates #######################
get_u_init <- function(elstH, N, K, d, lambda=1)
{
                                        # function to get initial estimate of the coordinates
                                        # want to use generalised MDS as in sarkar and moore
                                        # elstH is observed graph
                                        # N is number of nodes
                                        # K is max edge order
                                        # return (N x 2) matrix of latent coordinates
                                        # lambda is the scaling constant for the hyperedges

    ## convert hypergraph into a weighted matrix
    adj = matrix(0, nrow=N, ncol=N)
    for (i in 1:(N-1)){
        for (j in (i+1):N){
            idin = unlist(lapply( elstH, function(x){ sum(c(i,j) %in% x ) ==2 } ))
            if ( any( idin > 0 ) == TRUE ){
                adj[i,j] = 1
                adj[j,i] = adj[i,j]
            }
        }
    }
    
    ## convert to graph and calculate distances 
    grph = graph_from_adjacency_matrix( adj, mode="undirected", weighted=TRUE )
    dist = distances(grph, weights=E(grph)$weight)

    ## add in correction for disconnected nodes
    cc = components(grph)
    if (cc$no > 1){
        for (comp in cc$no:2){ # loop over smaller components and randmly connect       
            id_cc = which( cc$membership == comp )
            id_ncc = which( cc$membership != comp )
            ## randomly connect components
            n_cc = id_cc[sample( 1:length(id_cc), 1 )]
            n_ncc = id_ncc[sample( 1:length(id_ncc), 1 )]
            adj[n_cc,n_ncc] = mean(dist[ dist!=Inf] ) 
            adj[n_ncc,n_cc] = adj[n_cc,n_ncc] 
        }
	## update graph and distances
   	 grph = graph_from_adjacency_matrix( adj, mode="undirected", weighted=TRUE )
        dist = distances(grph, weights=E(grph)$weight)
    }
   
    ## ensure that there are no Inf distances (should be impossible!)
    if (sum(dist==Inf)>0){
        dist[dist==Inf] = max(dist[dist!=Inf]) + .2
    }
        
    ## get points from GMDS (assume points in R2)
    pts = cmdscale(dist, k=d) + mvrnorm(N, c(0,0), 0.01*diag(d))

    return(list(coords=pts, distances=dist))
}

######################### ABC init for mu and sig  #######################
get_mu_and_sig_init <- function(elstH, N, K, prior, r, phi0, phi1, idfixed, d, tol=10, nSmp = 20)
{
                                        # function to initialise mean and covariance of latent coords
                                        # use and ABC scheme to do this
                                        # elstH is observed graph
                                        # N is number of nodes
                                        # K is max e order
                                        # prior has hyperparameter values for mu and sigma
                                        # idfixed is the index of the fixed coordinates
                                        # tol is distance willing to accept, note this is discrete!

    ## calculate T for observed hypergraph
    lenH = unlist( lapply( elstH, length ) )
    nprs = sum(lenH==2)
    if (K > 2) {
        nhyp = sum(lenH==3)
    } else {
        nhyp = 0
    }

    ## calculate adj for observed graph
    adj = matrix(0, nrow=N, ncol=N)
    lenH = unlist( lapply(elstH, length) )
    prlst = elstH[lenH==2]
    nprs = length(prlst)
    if (nprs > 0) {
        for (i in 1:nprs){
            adj[prlst[[i]][1], prlst[[i]][2]] = 1
            adj[prlst[[i]][2], prlst[[i]][1]] = 1 # symmetric
        }
    }
    
    htmp = graph_from_adjacency_matrix( adj, mode="undirected" )
    ntri = length(triangles(htmp))/3
    TH = c( nprs, ntri, nhyp ) # summary statistics for  observed hypergraph

    muSmp = array(NA, c(nSmp, d))
    covSmp = array(NA, c(nSmp, d, d))
    uSmp = array(NA, c(N,d) )
    if (d == 2){
    uSmp[idfixed[1],] = c(-.5,0) # anchor points
    uSmp[idfixed[2],] = c(.5,0)
    } else if (d == 3){
        uSmp[idfixed[1],] = c(-.5,0,0) # anchor points
        uSmp[idfixed[2],] = c(.5,0,0)
    }
    s = 1
    while (s <= nSmp){

        ## propose mu and sigma
        muStr = mvrnorm(1, mu=prior$m, Sigma=prior$sigmu)
        sigStr = riwish( prior$nu, prior$Phi )

        ## sample nodes 
        uSmp[-idfixed,] = mvrnorm( N-length(idfixed), mu=muStr, Sigma=sigStr )

        ## get induced graph and perturb
        elstG = get_graph_from_cech(uSmp, r)
        elstG = addnoise(elstG, N, K, phi0, phi1)
        
        ## evaluate statistic T
        lenG = unlist( lapply( elstG, length ) )
        nprs = sum(lenG==2)
        if (K > 2){
            nhyp = sum(lenG==3)
        } else {
            nhyp = 0
        }

        ## get graph g
        adj = matrix(0, nrow=N, ncol=N)
        prlst = elstG[ lenG==2 ]
        nprs = length(prlst)
        if (nprs >= 1){
            for (i in 1:nprs){
                adj[prlst[[i]][1], prlst[[i]][2]] = 1
                adj[prlst[[i]][2], prlst[[i]][1]] = 1 # symmetric
            }
        }

        gtmp = graph_from_adjacency_matrix( adj, mode="undirected" )
        ntri = length(triangles(gtmp))/3
        TG = c( nprs, ntri, nhyp )

        ## accept/reject sample
        gdst = sqrt(sum(TH-TG)^2) # distance betwee summary statistics
        print( TH )
        print( TG )
        if ( gdst < tol ){
            print( paste("s = ", s ) )
            muSmp[s,] = muStr
            covSmp[s,,] = sigStr
            s = s+1
        } else {
            print( "reject proposal" )
        }
    }
    return( list( mus = muSmp, covs = covSmp ) )
}

########################### initialise radii #######################

get_r_init <- function(us, elstH, N, K, rInc){
                                        # function to get initial r
                                        # us is latent init
                                        # elstH is obs. graph
                                        # N is number of nodes
                                        # K is max e order

    ## # get rk init ###
    rInit = rep(0, K-1) # don't have r1
    lenH = unlist( lapply( elstH, length ) )
    for (k in 2:K)
    {
        ## find edges order k
        ek = elstH[ lenH==k ]
        nk = length(ek)
        ## apply miniball to each latent coordinate
        rstrk = rep(0, nk)
        if (nk > 0){
            for (i in 1:nk)
            {
                rstrk[i] = miniball_rsqr( t(us[ ek[[i]], ]), 2, k )
            }
            ## initialise at max rstr
            rInit[k-1] = quantile( sqrt( rstrk ), .4 )
        }
    }

    idzero = which(rInit==0)
    for (id in idzero){
        if ( id == 1 ){
            rInit[1]=0.1
        } else {
            rInit[id] = rInit[id-1] + .1
        }
    }

    ## sort to be ascending only if rInc is TRUE, otherwise rInc just needs to be positive
    if (rInc==TRUE){rInit = sort( rInit + runif( K-1, min=0, max=.01 )) }
    return( rInit )
}


############################### prior for sigma ##############################
get_sig_prior <- function( us )
{
                                        # function to get an estimate of prior covariance
                                        # does this using PCA on the bookstein coords
                                        # are the coordinates

    pca = prcomp( us ) # apply pca (need to do this again because of Bookstein)
    eigs = pca$sdev^2 # get eigenvalues
    cov = pca$rotation %*% diag(eigs) %*% t(pca$rotation)

    return( cov )
}

## script to simulate from multiple models and compare the resulting graphs
## will simulate from
## 1) our hypergraph model
## 2) model of brendan murphy and ng
## 3) beta model of feinberg



############################## LSH function ########################


################# clustering model ########################

clustmod <- function(a, phi, pi, tau, N, M){
                                        # function to simulate from the model of Ng and Murphy
                                        # see 'model-based clustering for random hypergraphs'
                                        # a is vector length K, with final entry = 1
                                        # phi is matrix size NxG
                                        # pi is a priori prob of edge topic label
                                        # tau is a priori prob od edge size label
                                        # N is number of nodes
                                        # M is number of edges

    elist = list() # stores edges as a list
    
    G = length( pi ) # number of topic classes
    K = length( tau ) # number of size classes
    a[K] = 1 # just to be sure constraint imposed
    
    xmat = matrix(NA, nrow=N, ncol=M)
    z1 = matrix(0, nrow=M, ncol=G) # topic label index for each edge
    z2 = matrix(0, nrow=M, ncol=K) # size label index for each edge
    for (m in 1:M) # loop over edges
    {
        ## sample topic label
        gTmp = sample( 1:G, 1, prob=pi )
        z1[m,gTmp] = 1 # record topic label
        
        ## sample size label
        kTmp = sample( 1:K, 1, prob=tau )
        z2[m,kTmp] = 1 # record size label

        ## for each i, sample whether or not it is in the edge
        ## now a[kTmp] * phi[ , gTmp] gives the prob of each node being in this edge
        for (n in 1:N)
        {
            xmat[n,m] = rbinom(1, 1, prob = a[kTmp] * phi[n,gTmp] )
        }

        ## find 1s for this edge and store in a list
        elist[[m]] = which(xmat[,m] != 0 )
    }
    
    return( list(X=xmat, Z1=z1, Z2=z2, elist=elist) )
}

################# beta model #############################

betamod <- function(beta, Kmax, unif = FALSE){
                                        # function to simulate from stasi et al
                                        # see \beta models for random hypergraphs with a given degree sequence
                                        # beta is the nodal parameter, vector length N
                                        # N is number of nodes
                                        # Kmax is maximum hyperedge order
                                        # unif indicates is model is k uniform
                                        # if unif = TRUE then Kmax denotes the order of all edges

    N = length(beta )
    
    ## output will be a list of edges (could convert into matrix later)
    edges = list()

    if (unif == FALSE){
        ## generate general hypergraph (no restriction on edge sizes) ##
        
        e = 1 # controls index for placement of edgelist
        for (k in 2:Kmax)
        {
            combs = combn(1:N, k) # all potential edges
            ncomb = dim(combs)[2] # number of edges
            for (i in 1:ncomb)
            {
                bTmp = beta[combs[,i]]
                pTmp = 1 / ( 1 + exp(-sum(bTmp)) )
                xTmp = rbinom( 1, 1, pTmp ) # flip for edge
                if (xTmp == 1)
                {
                    edges[[e]] = combs[,i] # add edge to edgelist
                    e = e + 1 # increment index
                }
            }
        }

    } else if (unif == TRUE) {
        ## generate k uniform hypergraph (with Kmax as edge size) ##

        e = 1 # placement for edgelist
        combs = combn(1:N, Kmax) # all potential edges
        ncomb = dim(combs)[2]
        for (i in 1:ncomb)
        {
            bTmp = beta[combs[,i]]
            pTmp = 1 / ( 1 + exp(-sum(bTmp)) )
            xTmp = rbinom( 1, 1, pTmp ) # flip for edge
            if (xTmp == 1)
            {
                edges[[e]] = combs[,i] # add edge to edgelist
                e = e + 1 # increment index
            }
        }
    }
    
    return( edges )
    
}

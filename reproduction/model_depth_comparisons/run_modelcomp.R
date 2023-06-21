library(praxi)

### script to run depth model comparisons


## call script with model functions
#source("compare_models.R") # script also contains testing for parameter values (commented in the region)

## define cases for each model
N = 50 
emax = 4

## BETA MODEL ##
prms_beta = list() # will be list of length 2
prms_beta[[1]] = list( beta=rep(-1.4, N), N=N )
prms_beta[[2]] = list( beta=seq(from=-0.5, to=-2, length=N), N=N )

## MURPHY MODEL ##
prms_murphy = list() # will be of length 4
prms_murphy[[1]] = list(a=c(1), phi = matrix( rep(.075, N), ncol=1 ), pi=c(1), tau=c(1), M=N^2, N=N)

typeind = sample( c(0,1,2), N, prob=rep(1/3,3), replace=TRUE ) 
phi = matrix( 0, nrow=N, ncol=3 )
phi[typeind==0,1] = .25
phi[typeind==1,2] = .25
phi[typeind==2,3] = .25
prms_murphy[[2]] = list(a=c(1), phi = phi, pi=c(1/3, 1/3, 1/3), tau=c(1), M=N^2, N=N)

prms_murphy[[3]] = list(a=c(.2, .5, 1), phi = matrix( .15, nrow=N, ncol=1 ), pi=c(1), tau=c(1/3, 1/3, 1/3), M=N^2, N=N)

typeind = sample( c(0,1,2), N, prob=rep(1/3,3), replace=TRUE ) # 0 = A, 1 = B, 2 = C
phi = matrix( 0, nrow=N, ncol=2 )
phi[typeind==0,1] = .3
phi[typeind==1,2] = .3
phi[typeind==2,] = rep(.2,2)
prms_murphy[[4]] = list(a=c(.4,1), phi = phi, pi=c(1/2,1/2), tau=c(1/3,1/3,1/3), M=N^2, N=N)

## LATENT SPACE MODEL ##
prms_lsm = list()
prms_lsm[[1]] = list( r=c(.18, .3, .35), phi0= c(.01, .01, .01), phi1=c(.01, .01, .01), m=c(0,0), sigm=.25*matrix( c(1, .9,.9, 1), nrow=2), N=N)
prms_lsm[[2]] = list( r=c(.18, .3, .35), phi0= c(.01, .01, .01), phi1=c(.01, .01, .01), m=c(0,0), sigm=.25*diag(2), N=N)
prms_lsm[[3]] = list( r=c(.2, .3, .35), phi0= c(.01, .01, .01), phi1= c(.01, .5, .01), m=c(0,0), sigm=.25*diag(2), N=N)
prms_lsm[[4]] = list( r=c(.1, .35, .4),  phi0= c(.01, .01, .01), phi1=c(.01, .01, .01), m=c(0,0), sigm=.25*diag(2), N=N)
prms_lsm[[5]] = list( r=c(.18, .3, .35),  phi0= c(.01, .01, .01), phi1=c(.01, .01, .01), m=c(0,0,0), sigm=.25*diag(3), N=N)

## ####### RUN SIMS ###################

## set up storage, seed
set.seed(123)
nReps = 100000 # number of graph simulations
## beta output
beta_out1 = list()
beta_out2 = list()
## murphy output
murphy_out1 = list()
murphy_out2 = list()
murphy_out3 = list()
murphy_out4 = list()
## lsm output
lsm_out1 = list()
lsm_out2 = list()
lsm_out3 = list()
lsm_out4 = list()
lsm_out5 = list()
## ksbm output
ksbm_out1 = list()
ksbm_out2 = list()
ksbm_out3 = list()
ksbm_out4 = list()
ksbm_out5 = list()


t = Sys.time()
## loop over repititions
for (rep in 1:nReps)
{

    print(rep)
    
    ## BETA MODEL ##
    
    ## loop over cases
    ncs_beta = length( prms_beta )
    for (c in 1:ncs_beta)
    {
        
        ## simulate graph
        out = betamod( prms_beta[[c]]$beta, emax )
        
        ## store
        if (c == 1)
        {
            beta_out1[[rep]] = getsummary(out, N)
        } else if (c ==2 ) {
            beta_out2[[rep]] = getsummary(out, N)
        }
    }
    print( "finished beta" )
    
    ## MURPHY MODEL ##

    ## loop over cases
    ncs_murphy = length( prms_murphy )
    for (c in 1:ncs_murphy)
    {
        
        ## simulate graph
        out = clustmod( prms_murphy[[c]]$a, prms_murphy[[c]]$phi, prms_murphy[[c]]$pi, prms_murphy[[c]]$tau, N, prms_murphy[[c]]$M)$elist
        
        ## store
        if (c == 1)
        {
            murphy_out1[[rep]] = getsummary(out, N)
        } else if (c ==2 ) {
            murphy_out2[[rep]] = getsummary(out, N)
        } else if (c==3) {
            murphy_out3[[rep]] = getsummary(out, N)
        } else if (c==4) {
            murphy_out4[[rep]] = getsummary(out, N)
        }
    }
    print( "finished murphy")
    
    ## LS MODEL ##

    ## loop over cases
    ncs_lsm = length( prms_lsm )
    for (c in 1:ncs_lsm)
    {
        
        ## simulate graph
        out = genHypergraph(N, emax, prms_lsm[[c]]$r, prms_lsm[[c]]$phi0, prms_lsm[[c]]$phi1, prms_lsm[[c]]$m, prms_lsm[[c]]$sigm, add_noise=TRUE, bkstn=FALSE)$elst
        
        ## store
        if (c == 1)
        {
            lsm_out1[[rep]] = getsummary(out, N)
        } else if (c ==2 ) {
            lsm_out2[[rep]] = getsummary(out, N)
        } else if (c==3) {
            lsm_out3[[rep]] = getsummary(out, N)
        } else if (c==4) {
            lsm_out4[[rep]] = getsummary(out, N)
        } else if (c==5) {
            lsm_out5[[rep]] = getsummary(out, N)
        }
    }
    print( "finished lsm" )

    if ( rep %% 50 == 0 ){print('#')}
    if ( rep %% 1000 == 0)
    {

        ## store output
        save(beta_out1, beta_out2, file="./output/beta_out.RData") # variables will be loaded in with the same names
        save(murphy_out1, murphy_out2, murphy_out3, murphy_out4, file="./output/murphy_out.RData")
        save(lsm_out1, lsm_out2, lsm_out3, lsm_out4, lsm_out5, file="./output/lsm_out.RData")
        
    }
}

t = Sys.time() - t

## store output
save(beta_out1, beta_out2, file="./output/beta_out.RData") # variables will be loaded in with the same names
save(murphy_out1, murphy_out2, murphy_out3, murphy_out4, file="./output/murphy_out.RData")
save(lsm_out1, lsm_out2, lsm_out3, lsm_out4, lsm_out5, file="./output/lsm_out.RData")

## store parameters
save(prms_beta, prms_murphy, prms_lsm, file="./output/prms.RData")

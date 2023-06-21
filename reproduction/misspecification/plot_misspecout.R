## summarise output for misspecification simulation study
library(ggplot2)
library(dplyr)
library('HyperG')
ibrary(hyperlatent)
#source('../mcmc/mcmc.R')
#source('../mcmc/getdat.R')
#source("../mcmc/graphsummary.R")
#source("../mcmc/getmisspecdat.R") # functions to gen misspec data
#source("../mcmc/misspecsummary.R") # function getting graph summaries

## extra functions
getavdd <- function(ddlst){

    nrep = length(ddlst)
    maxdd = max( unlist( lapply(ddlst, length) ) )
    avdd = rep(0, maxdd)
    for (i in 1:nrep){
        ddtmp = ddlst[[i]]
        ltmp = length(ddtmp)
        avdd[1:ltmp] = avdd[1:ltmp] + ddtmp
    }
    avdd = avdd/nrep

    return(avdd)
}
plotddn <- function(ddtrue, ddpost){
                                        # function to compare degree distributions
    
    ## ddest is estimate of dd
    xmax = max( max(ddtrue), max(ddpost) )
    
    par(mfrow=c(1,2), oma=c(0,0,2,0))
    par(mar=c(4,4,1,1))
    barplot(ddtrue, xlab="Degree", ylab="Probability",  names.arg=0:(length(ddtrue)-1), ylim=c(0, max(ddtrue)), cex.lab=1.25, cex.axis=1.2, cex.main=1.25 )
    qqplot( ddtrue, ddpost, xlab="'Truth'", ylab="Estimate", xlim = c(0, xmax), ylim=c(0,xmax), cex=1.2, pch=16 ,cex.lab=1.25, cex.axis=1.25, cex.main=1.25 )
    abline(0,1, lwd = 3, col='red3')   
}

#######################################################################
## cases are
## umiss: 1 none, 2 clust1, 3 clust2, 4 uniform
## rmiss: 1 none, 2 simplicial
## phimiss: 1 none
umiss=c(1,2,3,4)
rmiss=c(1,2)
phimiss=c(1)
cases = do.call( expand.grid, list(u=umiss, r=rmiss, phi=phimiss))
ncs = dim(cases)[1]

## loop over each case and check output
for (ctmp in 1:ncs){
    ucase = cases[ctmp,1]
    rcase = cases[ctmp,2]
    phicase = cases[ctmp,3]
    label = paste("prms_u",ucase,"r",rcase,"phi",phicase, sep="")
    ## read in data
    load( file=paste("./output/H_", label,  ".RData", sep="") ) # H
    ## read in output
    load( file=paste("./output/mcmcout_", label,  ".RData", sep="") ) # mcmcout
    strg = mcmcout$strg
    
    nIts = dim(strg$uTrc)[3]
    nBrn = ceiling(nIts/2)
    indx = seq(nBrn, nIts, by=10) 
     
    ## ############# INSPECT OUTPUT #########################
    ## plot output for sanity
    png(file=paste("./output/us_trc_", label,  ".png", sep=""), width=400, height=400)
    par(mfrow=c(1,1))
    plot_mcmc_output_us_only(strg, indx, H$elst, inclH=FALSE)
    dev.off()

    png(file=paste("./output/nonus_trc_", label,  ".png", sep=""), width=1200, height=400)
    plot_mcmc_output_non_us_trcs(strg, indx, H$elst)
    dev.off()

    ## ############# GET PREDICTIVES #########################
    load( file=paste("./output/pred_summary_", label,  ".RData", sep="") ) # preds
    load( file=paste("./output/true_summary_", label,  ".RData", sep="") ) # truesum
    
    ## ############# PLOT THE PREDICTIVES #########################

    ## have 1) degree distribution, 2) motifs: ntri, h1, h2, h3, ijk, 3) clustering
    png( file=paste("./output/misspec_summ_", label,  ".png", sep=""), width=1000, height=175)
    par(mfrow=c(1,7))
    par(oma = c(4,4,1,1))
    par(mar = c(4,4,1,1))
    ## degree distribution
    dd_smp_av = getavdd(preds$DD)
    qqplot(dd_smp_av, truesum$DD, xlab="Estimate", ylab="Truth", cex.lab=1.2, cex.axis=1.4)
    abline(0,1, lwd=3, col='red')
    ## motif counts
    boxplot(preds$ntri, xlab="m3", ylab="count", boxlwd=3, cex.lab=1.4, cex.axis=1.4) 
    abline(h=truesum$ntri, col='red', lwd=3)
    boxplot(unlist(preds$h1), xlab="h1", ylab="count", boxlwd=3, cex.lab=1.4, cex.axis=1.4) 
    abline(h=truesum$h1, col='red', lwd=3)
    boxplot(unlist(preds$h2), xlab="h2", ylab="count", boxlwd=3, cex.lab=1.4, cex.axis=1.4) 
    abline(h=truesum$h2, col='red', lwd=3)
    boxplot(unlist(preds$h3), xlab="h3", ylab="count" , boxlwd=3, cex.lab=1.4, cex.axis=1.4)
    abline(h=truesum$h3, col='red', lwd=3)
    boxplot(preds$ijk, xlab="ijk", ylab="count", boxlwd=3, cex.lab=1.4, cex.axis=1.4) 
    abline(h=truesum$ijk, col='red', lwd=3)
    ## clustering
    boxplot(preds$nclst, xlab="nclst", ylab="count", boxlwd=3, cex.lab=1.4, cex.axis=1.4) 
    abline(h=truesum$nclst, col='red', lwd=3)
    dev.off()

}

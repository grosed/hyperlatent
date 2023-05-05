library(reshape2)
library(ggplot2)
library('HyperG') 
source("../mcmc/graphsummary.R")
source('../mcmc/mcmc.R')
source('../mcmc/getdat.R')

########################################################################
################# extra functions #####################################
########################################################################

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

convert_list_to_boxplot_df <- function(l){
    ## need two vectors 1) contains probabilities, 2) contains degree
    probs = c()
    degs = c()
    for (i in 1:length(l)){
        probs = append(probs, l[[i]])
        degs = append(degs, 0:(length(l[[i]])-1))
    }
    df = data.frame(probs=probs, degs=as.factor(degs))
    return(df)
}

convert_joint_list_to_boxplot_df <- function(lest, ltru){
    ## need two vectors 1) contains probabilities, 2) contains degree
    probs_est = c()
    degs_est = c()
    ## add on posterior estimates
    for (i in 1:length(lest)){
        probs_est = append(probs_est, lest[[i]])
        degs_est = append(degs_est, 0:(length(lest[[i]])-1))
    }
    probs_tru = c()
    degs_tru = c()    
    ## add on true 
    for (i in 1:length(ltru)){
        probs_tru = append(probs_tru, ltru[[i]])
        degs_tru = append(degs_tru, 0:(length(ltru[[i]])-1))
    }
    type = c(rep("post", length(probs_est)), rep("true", length(probs_tru)))
    df = data.frame(probs=c(probs_est,probs_tru), degs=as.factor(c(degs_est, degs_tru)), type=as.factor(type))
    return(df)
}

convert_list_to_upperlower_df <- function(bxplt_df){
    ## make this from the bxplot_df
    ## for each probability have lower, upper, median
    ds = sort(unique(bxplt_df$degs))
    uppr = rep(NA, length(ds))
    lwr = rep(NA, length(ds))
    medn = rep(NA, length(ds))
    for (i in 1:length(ds)){
        ptmp = df$probs[ df$degs == ds[i] ]
        uppr[i] = max(ptmp)
        lwr[i] = min(ptmp)
        medn[i] = median(ptmp)
    }
    df = data.frame(degs=as.factor(ds), uppr=uppr, lwr=lwr, medn=medn)
    return(df)
}

########################################################################
########################################################################

## load in data
load( file="./output/postfit.RData") # strg
load( file="./output/post_time.RData") # time
load( file="./output/H.RData") # H

## plot the data 
png(file="./plots/H_S_6_1_data.png", width=400, height=400)
plotHypergraph(H$us, H$elst, K)
dev.off()

## check the mcmc fit
nIts = dim(strg$ll)[2]
indx = seq(75000, nIts, by=20)
plot_mcmc_output(strg, indx, H$elst) 
dev.off() 

png(file="./plots/H_S_6_1_us_fit.png", width=800, height=400)
plot_mcmc_output_us_only(strg, indx, H$elst)
dev.off()

png(file="./plots/H_S_6_1_all_but_us_fit.png", width=2500, height=1250, res=300)
plot_mcmc_output_non_us_trcs(strg, indx, H$elst)
dev.off()

## true parameters
prmsTru = list(us=H$us, r = c(.35, .45), phi0 = rep(.0001, 2), phi1 = rep(.0001, 2), N=50, K=3)

prmsEst = list(us=apply(strg$uTrc[,,indx], c(1,2), mean), r = apply(strg$rTrc[,indx], 1, mean), phi0 = apply(strg$phi0Trc[,indx], 1, mean), phi1 = apply(strg$phi1Trc[,indx], 1, mean), mu=apply(strg$muTrc[,indx], 1, mean), sigma=apply(strg$sigTrc[,,indx], c(1,2), mean))

## compare degree distributions
load(file="./output/true_dd.RData") # true_dd
load(file="./output/true_dd_smp.RData") # true_dd_smp
load(file="./output/true_dd_smp_av.RData") # true_dd_smp_av
load(file="./output/postout_dd.RData") # post_dd
load(file="./output/post_dd_av.RData") # post_dd_av

nrep = length(post_dd$full)

xwd = 6
ywd = 3

#########################################################################
## make plots
#########################################################################

## make 'old' boxplot connections
pdf(file="./plots/S_6_1_truevspost_old_minmax.pdf", width=xwd, height=ywd)
df = convert_list_to_boxplot_df(post_dd$old)
df = convert_list_to_upperlower_df(df)
ggplot(data = df, mapping = aes(x = degs, y = medn)) + geom_pointrange(mapping = aes(ymin = lwr, ymax = uppr)) + geom_point( data=data.frame(x=0:(length(true_dd$full)-1)+1, y=true_dd$full), aes(x=x,y=y), color='green4', shape=17, size=3) + theme_classic() + labs(x="degree", y="probability") + theme(axis.text.x = element_text(angle=90)) + scale_x_discrete(breaks=levels(df$degs)[c(T, rep(F,9))]) 
dev.off()

pdf(file="./plots/S_6_1_truevspost_oldk2_minmax.pdf", width=xwd, height=ywd)
df = convert_list_to_boxplot_df(post_dd$oldk2)
df = convert_list_to_upperlower_df(df)
ggplot(data = df, mapping = aes(x = degs, y = medn)) + geom_pointrange(mapping = aes(ymin = lwr, ymax = uppr)) + geom_point( data=data.frame(x=0:(length(true_dd$fullk2)-1)+1, y=true_dd$fullk2), aes(x=x,y=y), color='green4', shape=17, size=3) + theme_classic() + labs(x="degree", y="probability") + theme(axis.text.x = element_text(angle=90)) + scale_x_discrete(breaks=levels(df$degs)[c(T, rep(F,4))]) 
dev.off()

pdf(file="./plots/S_6_1_truevspost_oldk3_minmax.pdf", width=xwd, height=ywd)
df = convert_list_to_boxplot_df(post_dd$oldk3)
df = convert_list_to_upperlower_df(df)
ggplot(data = df, mapping = aes(x = degs, y = medn)) + geom_pointrange(mapping = aes(ymin = lwr, ymax = uppr)) + geom_point( data=data.frame(x=0:(length(true_dd$fullk3)-1)+1, y=true_dd$fullk3), aes(x=x,y=y), color='green4', shape=17, size=3) + theme_classic() + labs(x="degree", y="probability") + theme(axis.text.x = element_text(angle=90)) + scale_x_discrete(breaks=levels(df$degs)[c(T, rep(F,9))]) 
dev.off()

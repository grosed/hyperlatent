source('../mcmc/getdat.R')
source('../mcmc/mcmc.R')
source("../mcmc/graphsummary.R") 
source("./data/subsample.R")
library('HyperG')
library('ggplot2')
library('gridExtra')

#################################### additional functions ################################
plot_hypergraph_w_hyperg <- function(elst, us){
    N = dim(us)[1]
    h = hypergraph_from_edgelist(elst, v=1:N)
    plot(h, layout=us)
}

summ_pred_dist_tru <- function(est_df, true_df){
                                        # create list to summarise 'distance' between true and predictive distn
    
    ## output as a list
    out = list()
    
    ## get list of nodes and types to loop through
    node_ids = levels(est_df$node)
    type_ids = levels(est_df$type)

    ## loop over types (different motifs)
    for (ty in type_ids){

        ltmp = list()

        for (n in node_ids){
            ## find distances
            dtmp = subset(est_df, node == n & type == ty)$count - subset(true_df, node == n & type == ty )$count
            ltmp[[n]] = abs(dtmp) # take absolute distance
        }
        dtmp = unlist(ltmp) # create 'master' list of all predictive nodes
        
        max_dist = max(dtmp) # largest discrepancy
        min_dist = min(dtmp) # smallest discrepancy (can be < 0)
        dist_vec = seq(min_dist, max_dist, by=1)
        
        ptmp = rep(NA, length(dist_vec))
        for (p in 1:length(dist_vec)){ ptmp[p] = sum( dtmp== dist_vec[p] ) / length(dtmp) }

        out[[ty]]$ps = ptmp
        out[[ty]]$ds = dist_vec
    }

    return( out )
}

#########################################################################
################################ DBLP ###################################
#########################################################################

## read in data to fit
load("./data/dblp_rw_to_fit.RData")
elst = rw_sub_elst
N = length(unique(unlist(elst)))
K = max( unlist( lapply( elst, length ) ) )

## plots of data using HyperG hypergraph function
dblp = hypergraph_from_edgelist(elst, v=1:N)
 pdf("./plots/dblp_hypergraph.pdf", width=10, height=10)
plot.hypergraph( dblp, layout = layout.graphopt, vertex.size=7, mark.col=rainbow(dim(dblp$M)[1], alpha=.2))
dev.off()

## plot random walk subset
pdf("./plots/dblp_hypergraph_rw_subsample.pdf", width=10, height=10)
set.seed(633)
rw_ndsmp = random_walk_subsample(elst, 5, c_ret = 0.05, nmin = 30)
dblp_subsmp = induced_hypergraph(dblp, rw_ndsmp )
plot.hypergraph( remove.isolates(dblp_subsmp), layout = layout.graphopt, vertex.size=7, mark.col=rainbow(dim(dblp$M)[1], alpha=.2))
dev.off()

## read in 'full' data (ie with additional nodes)
load("./data/dblp_rw_full.RData")
elst_full = rw_full_elst

## locate the 'new nodes'
nd_new = setdiff( unique(unlist(elst_full)), unique(unlist(elst)) )
Nstr = length(nd_new)

## find the summary of the 'true' additional nodes
summ_true = getsummary_for_node_subset(elst_full, nd_new)
save(summ_true, file="./output/dblp_summ_true_new_nodes.RData")

hypmotif_true = get_hypmotif_count(elst_full)
save(hypmotif_true, file="./output/dblp_hypmotif_true.RData")

## calculate predictives 
set.seed(57234)
for (ID in c("simp", "nonsimp")){
    
    load( file=paste("./output/dblp_rw_", ID, "_fit.RData", sep="")) 

    ## check the output of the fitted model
    nIts = dim(out$strg$ll)[2]
    nBrn = ceiling( 3 * nIts / 4 )
    ids = seq(nBrn, nIts, 20) 
    nRep = length(ids)

    ## inspect hypergraph from posterior mean
    rpst = apply(out$strg$rTrc[,ids],1,mean)
    upst = apply( out$strg$uTrc[,,ids], c(1,2), mean )
    phi0pst = apply(out$strg$phi0Trc[,ids],1,mean)
    phi1pst = apply(out$strg$phi1Trc[,ids],1,mean)
    
    png(paste("./plots/dblp_", ID, "_fit_hyp_upst.png", sep=""), width=5000, height=500)
    par(mfrow=c(1,1))
    plot_hypergraph_w_hyperg(elst, upst)
    dev.off() 

    png(paste("./plots/dblp_", ID, "_fit_nonus_trc.png", sep=""), width=1000, height=400)
    plot_mcmc_output_non_us_trcs(out$strg, ids, elst, xlims = ?, ylims = ?)
    dev.off()
    
    png(paste("./plots/dblp_", ID, "_fit_us_trc.png", sep=""), width=600, height=600)
    plot_mcmc_output_us_only(out$strg, ids, elst, inclH=FALSE)
    dev.off()

    png(paste("./plots/dblp_", ID, "_fit_us_trc_zoomed.png", sep=""), width=600, height=600)
    plot_mcmc_output_us_only(out$strg, ids, elst, inclH=FALSE, xlims= c(-1.5,2.5), ylims =c(-4,3) )
    dev.off()   
    
    ## summarise posterior predictive for new nodes, and for hypergraph motif counts    
    count_types = c("degk2", "degk3", "degk4", "tris", "n1", "n2", "n3")
    count = rep(NA, Nstr*length(count_types)*nRep)
    type = rep(rep(count_types, each=Nstr), nRep)
    node = rep(rep(1:Nstr, length(count_types)), nRep)

    ## posterior predictives for counting number of non-simplicial motifs
    krng = 3:K
    count_ns_vs_s = rep(NA, Nstr*nRep*2*length(krng))
    ns_vs_s_type = rep(rep(c("nonsimp", "simp"), each=Nstr*length(krng), nRep))
    ns_vs_s_node = rep(rep(rep(1:Nstr, length(krng)),2), nRep)
    ns_vs_s_k = rep(rep(rep(krng, each=Nstr), 2), nRep)    
    
    hypmotif_types = c("n1", "n2", "n3")
    hypmotif_count = rep(NA, nRep*length(hypmotif_types))

    set.seed(7845)
    for( rp in 1:nRep ){
        
        ## sample hypergraph, gives eOld, eNew, eBtwn
        prms = list(N = N, K = K, us = out$strg$uTrc[,,ids[rp]], mu = out$strg$muTrc[,ids[rp]], sigma = out$strg$sigTrc[,,ids[rp]], r = out$strg$rTrc[,ids[rp]], phi0 = out$strg$phi0Trc[,ids[rp]], phi1 = out$strg$phi1Trc[,ids[rp]])
        htmp = genHyperG_with_NewConnection(Nstr, prms, periph=TRUE, addnoise=TRUE)

        ## combine list into full list
        elst_tmp = c(htmp$eOld, htmp$eNew, htmp$eBtwn) # the full elist

        ## summary of connections for new nodes
        summ_tmp = getsummary_for_node_subset(elst_tmp, nd_new)
        count[((rp-1)*length(count_types)*Nstr + 1):(rp*length(count_types)*Nstr)] = unlist(summ_tmp)

        ## count number of simplicial and nonsimplicial motifs each new node belongs to
        ns_count = count_orderk_nonsimp_motifs(elst_tmp, krng, nd_new)
        count_ns_vs_s[((rp-1)*(Nstr*length(krng)*2) + 1):( rp*(Nstr*length(krng)*2) )] = c(c(t(ns_count$count_ns)), c(t(ns_count$count_s))) 
        
        ## sub-hypergraph count (for full hypergraph)
        hypmotif_count[((rp-1)*length(hypmotif_types) + 1):(rp*length(hypmotif_types)) ] = unlist(get_hypmotif_count(elst_tmp))
        
        if( rp %% 100 == 0 ){print(paste("Finished rep ", rp))}
    }

    ## posterior predictive summary for new nodes
    summ = data.frame(count=count, node=as.factor(node), type=as.factor(type))
    save(summ, file=paste("./output/dblp_", ID, "_summ_pred_new_nodes.RData", sep=""))

    ## posterior predictive membership of nonsimplicial motifs for new nodes
    summ_ns = data.frame(count=count_ns_vs_s, node=as.factor(ns_vs_s_node), type=as.factor(ns_vs_s_type), k=as.factor(ns_vs_s_k))
    save(summ_ns, file=paste("./output/dblp_", ID, "_summ_nsnmotif_pred_new_nodes.RData", sep=""))
    
    ## posterior pred motif counts full hypergraph
    pred_motif = data.frame(count=hypmotif_count, type=as.factor(rep(hypmotif_types, nRep)))
    save(pred_motif, file=paste("./output/dblp_", ID, "_motif_predictive.RData", sep=""))

    ## summarise predictive 'distance from the truth'
    summ_tru = data.frame(count=unlist(summ_true), node=factor(rep(1:Nstr, length(count_types))), type=factor(rep(count_types, each=Nstr)))
    diff_dblp = summ_pred_dist_tru(summ, summ_tru)
    save(diff_dblp, file=paste("./output/dblp_", ID, "_diff_pred_from_true.RData", sep=""))
        
    ## plot with ggplot
    cols = c("blue", "purple", "cyan4", "forestgreen", "black", "orchid", "darkorange")

    pdf(paste("./plots/dblp_", ID, "_pred_summary.pdf", sep=""), width=10, height=6)

    p1 = ggplot() + geom_boxplot(data=subset(summ, node %in% 1:5), aes(x=node, y=count, color=type), fill='grey80') + scale_color_manual(values=cols) + theme_minimal() + geom_boxplot(data=subset(summ_tru, node %in% 1:5), aes(x=node, y=count, fill=type), color='red', lwd=1.5)+ guides(colour = guide_legend(override.aes = list(color = cols), lwd=1, fill='grey80'))

    p2 = ggplot() + geom_boxplot(data=subset(summ, node %in% 6:10), aes(x=node, y=count, color=type), fill='grey80') + scale_color_manual(values=cols) + theme_minimal() + geom_boxplot(data=subset(summ_tru, node %in% 6:10), aes(x=node, y=count, fill=type), color='red', lwd=1.5)+ theme(legend.position="none")

    grid.arrange(p1, p2, nrow=2)

    dev.off()
}

#########################################################################
######### Make histogram to report predictive accuracy

ID = "nonsimp"
load(file=paste("./output/dblp_", ID, "_diff_pred_from_true.RData", sep="")) # diff_dblp

pdf(file=paste("./plots/dblp_",ID, "_pred_bp_degree.pdf", sep=""), width=6, height=2)
par(mfrow=c(1,3))
barplot(diff_dblp$degk2$ps, names.arg=diff_dblp$degk2$ds, xlab="D", main="k=2", ylab="Proportion")
barplot(diff_dblp$degk3$ps, names.arg=diff_dblp$degk3$ds, xlab="D", main="k=3", ylab="Proportion")
barplot(diff_dblp$degk4$ps, names.arg=diff_dblp$degk4$ds, xlab="D", main="k=4", ylab="Proportion")
dev.off()

pdf(file=paste("./plots/dblp_",ID, "_pred_bp_motifs.pdf", sep=""), width=8, height=2)
par(mfrow=c(1,4))
barplot(diff_dblp$n1$ps, names.arg=diff_dblp$n1$ds, xlab="D", main="h1", ylab="Proportion")
barplot(diff_dblp$n2$ps, names.arg=diff_dblp$n2$ds, xlab="D", main="h2", ylab="Proportion")
barplot(diff_dblp$n3$ps, names.arg=diff_dblp$n3$ds, xlab="D", main="h3", ylab="Proportion")
barplot(diff_dblp$tris$ps, names.arg=diff_dblp$tris$ds, xlab="count", main="m3", ylab="Proportion")
dev.off()

#########################################################################
######### COMPARE SIMPLICIAL AND NON SIMPLICIAL ON SUBSET

## subset the data
nd_cntr = 27
id_ndin = unlist(lapply(elst, function(x){ nd_cntr %in% x }))
elst_nd_sub = elst[id_ndin] 
nd_sub = sort(unique(unlist(elst_nd_sub)))
h_sub = induced_hypergraph(dblp, nd_sub)
elst_sub = unique(hypergraph_as_edgelist(h_sub))
h_sub = hypergraph_from_edgelist(elst_sub)
plot(h_sub)

png(file="./plots/dblp_subgraph_comp.png", width=900, height=300)
par(mfrow=c(1,3)) 
plot(h_sub) 
load( file="./output/dblp_rw_simp_fit.RData")
simp_out = out
load( file="./output/dblp_rw_nonsimp_fit.RData")
nonsimp_out = out
nIts = dim(simp_out$strg$ll)[2]
nBrn = ceiling( 3 * nIts / 4 )
ids = seq(nBrn, nIts, 20) 
plot_mcmc_output_us_only(simp_out$strg, ids, elst, nd_sub, inclH=FALSE, jitr=.01)
plot_mcmc_output_us_only(nonsimp_out$strg, ids, elst, nd_sub, inclH=FALSE, jitr=.01)
dev.off()

#########################################################################
############################# GROCERY ###################################
#########################################################################

## read in data to fit
load("./data/gro_rndm_to_fit.RData")
elst = rndm_sub_elst
N = length(unique(unlist(elst)))
K = max( unlist( lapply( elst, length ) ) )

## plots of data using HyperG as hypergraph
gro = hypergraph_from_edgelist(elst)
pdf("./plots/gro_hypergraph.pdf", width=10, height=10)
plot( gro) 
dev.off()

## plot random walk subset
pdf("./plots/gro_hypergraph_rw_subsample.pdf", width=10, height=10)
set.seed(633)
rw_ndsmp = random_walk_subsample(elst, 5, c_ret = 0.05, nmin = 20)
gro_subsmp = induced_hypergraph(gro, rw_ndsmp )
plot.hypergraph( remove.isolates(gro_subsmp), layout = layout.graphopt, vertex.size=7, mark.col=rainbow(dim(gro$M)[1], alpha=.2))
dev.off()

## read in 'full' data (ie with additional nodes)
load("./data/gro_rndm_full.RData")
elst_full = rndm_full_elst

## locate the 'new nodes'
nd_new = setdiff( unique(unlist(elst_full)), unique(unlist(elst)) )
Nstr = length(nd_new)

## find the summary of the 'true' additional nodes
summ_true = getsummary_for_node_subset(elst_full, nd_new)
save(summ_true, file="./output/grocery_summ_true_new_nodes.RData")

hypmotif_true = get_hypmotif_count(elst_full)
save(hypmotif_true, file="./output/grocery_hypmotif_true.RData")

## calculate predictives 
set.seed(345)
for (ID in c("simp", "nonsimp")){ #, "rnoninc")){

    load( file=paste("./output/grocery_fit_", ID, "_rndm.RData", sep="" ))

    ## check the output of the fitted model
    nIts = dim(out$strg$ll)[2]
    nBrn = (3 * nIts) / 4 
    ids = seq(nBrn, nIts, 10) 
    nRep = length(ids)
    rpst = apply(out$strg$rTrc[,ids],1,mean)
    upst = apply( out$strg$uTrc[,,ids], c(1,2), mean )
    
    png(paste("./plots/gro_fit_", ID, "_hyp_upst.png", sep=""), width=500, height=500)
    uPst = apply( out$strg$uTrc[,,ids], c(1,2), mean )
    plot_hypergraph_w_hyperg(elst, upst)
    dev.off()

    png(paste("./plots/gro_fit_", ID, "_nonus_trc.png", sep=""), width=1000, height=400)
    plot_mcmc_output_non_us_trcs(out$strg, ids, elst) 
    dev.off()

    png(paste("./plots/gro_fit_", ID, "_us_trc.png", sep=""), width=600, height=600)
    plot_mcmc_output_us_only(out$strg, ids, elst, inclH=FALSE)
    dev.off()

    png(paste("./plots/gro_fit_", ID, "_us_trc_zoomed.png", sep=""), width=600, height=600)
    plot_mcmc_output_us_only(out$strg, ids, elst, inclH=FALSE, xlims= c(-5,5), ylims =c(-5,5) )
    dev.off()
    
    ## summarise posterior predictive for new nodes, and for hypergraph motif counts    
    count_types = c("degk2", "degk3", "degk4", "tris", "n1", "n2", "n3")
    count = rep(NA, Nstr*length(count_types)*nRep)
    type = rep(rep(count_types, each=Nstr), nRep)
    node = rep(rep(1:Nstr, length(count_types)), nRep)

    ## posterior predictives for counting number of non-simplicial motifs
    krng = 3:K
    count_ns_vs_s = rep(NA, Nstr*nRep*2*length(krng))
    ns_vs_s_type = rep(rep(c("nonsimp", "simp"), each=Nstr*length(krng), nRep))
    ns_vs_s_node = rep(rep(rep(1:Nstr, length(krng)),2), nRep)
    ns_vs_s_k = rep(rep(rep(krng, each=Nstr), 2), nRep)

    ## storage for storing motifs for full hypergraph
    hypmotif_types = c("n1", "n2", "n3")
    hypmotif_count = rep(NA, nRep*length(hypmotif_types))

    set.seed(234)
    for( rp in 1:nRep ){
        
        ## sample hypergraph, gives eOld, eNew, eBtwn
        prms = list(N = N, K = K, us = out$strg$uTrc[,,ids[rp]], mu = out$strg$muTrc[,ids[rp]], sigma = out$strg$sigTrc[,,ids[rp]], r = out$strg$rTrc[,ids[rp]], phi0 = out$strg$phi0Trc[,ids[rp]], phi1 = out$strg$phi1Trc[,ids[rp]])
        htmp = genHyperG_with_NewConnection(Nstr, prms, periph=TRUE, addnoise=TRUE)

        ## combine list into full list
        elst_tmp = c(htmp$eOld, htmp$eNew, htmp$eBtwn) # the full elist
     
        ## summary of connections for new nodes (degree and hypergraph motif counts)
        summ_tmp = getsummary_for_node_subset(elst_tmp, nd_new)
        count[((rp-1)*length(count_types)*Nstr + 1):(rp*length(count_types)*Nstr)] = unlist(summ_tmp)  
        ## count number of simplicial and nonsimplicial motifs each new node belongs to
        ns_count = count_orderk_nonsimp_motifs(elst_tmp, krng, nd_new)
        count_ns_vs_s[((rp-1)*(Nstr*length(krng)*2) + 1):( rp*(Nstr*length(krng)*2) )] = c(c(t(ns_count$count_ns)), c(t(ns_count$count_s))) 

        ## sub-hypergraph count for full hypergraph
        hypmotif_count[((rp-1)*length(hypmotif_types) + 1):(rp*length(hypmotif_types)) ] = unlist(get_hypmotif_count(elst_tmp))

        if( rp %% 100 == 0 ){print(paste("Finished rep ", rp))}
    }

    ## posterior predictive summary for new nodes
    summ = data.frame(count=count, node=as.factor(node), type=as.factor(type))
    save(summ, file=paste("./output/gro_", ID, "_summ_pred_new_nodes.RData", sep=""))

    ## posterior predictive membership of nonsimplicial motifs for new nodes
    summ_ns = data.frame(count=count_ns_vs_s, node=as.factor(ns_vs_s_node), type=as.factor(ns_vs_s_type), k=as.factor(ns_vs_s_k))
    save(summ_ns, file=paste("./output/gro_", ID, "_summ_nsnmotif_pred_new_nodes.RData", sep=""))

    ## posterior pred motif counts full hypergraph
    pred_motif = data.frame(count=hypmotif_count, type=as.factor(rep(hypmotif_types, nRep)))
    save(pred_motif, file=paste("./output/gro_", ID, "_motif_predictive.RData", sep=""))

    ## summarise predictive 'distance from the truth'    
    summ_tru = data.frame(count=unlist(summ_true), node=factor(rep(1:Nstr, length(count_types))), type=factor(rep(count_types, each=Nstr)))
    diff_gro = summ_pred_dist_tru(summ, summ_tru)
    save(diff_gro, file=paste("./output/gro_", ID, "_diff_pred_from_true.RData", sep=""))
    
    ## plot with ggplot
    cols = c("blue", "purple", "cyan4", "forestgreen", "black", "orchid", "darkorange")

    pdf(paste("./plots/gro_", ID, "_pred_summary.pdf", sep=""), width=15, height=4)
    ggplot() + geom_boxplot(data=summ, aes(x=node, y=count, color=type), fill='grey80') + scale_color_manual(values=cols) + theme_minimal() + geom_boxplot(data=summ_tru, aes(x=node, y=count, fill=type), color='red', lwd=1.5) + guides(colour = guide_legend(override.aes = list(color = cols), lwd=1, fill='grey80'))
    dev.off()
}

#########################################################################
######### Make histogram to report predictive accuracy

ID = "nonsimp"
load(file=paste("./output/gro_", ID, "_diff_pred_from_true.RData", sep=""))

pdf(file=paste("./plots/gro_",ID, "_pred_bp_degree.pdf", sep=""), width=6, height=2)
par(mfrow=c(1,3))
barplot(diff_gro$degk2$ps, names.arg=diff_gro$degk2$ds, xlab="D", main="k=2", ylab="Proportion")
barplot(diff_gro$degk3$ps, names.arg=diff_gro$degk3$ds, xlab="D", main="k=3", ylab="Proportion")
barplot(diff_gro$degk4$ps, names.arg=diff_gro$degk4$ds, xlab="D", main="k=4", ylab="Proportion")
dev.off()

pdf(file=paste("./plots/gro_",ID, "_pred_bp_motifs.pdf", sep=""), width=8, height=2)
par(mfrow=c(1,4))
barplot(diff_gro$n1$ps, names.arg=diff_gro$n1$ds, xlab="D", main="h1", ylab="Proportion")
barplot(diff_gro$n2$ps, names.arg=diff_gro$n2$ds, xlab="D", main="h2", ylab="Proportion")
barplot(diff_gro$n3$ps, names.arg=diff_gro$n3$ds, xlab="D", main="h3", ylab="Proportion")
barplot(diff_gro$tris$ps, names.arg=diff_gro$tris$ds, xlab="D", main="m3", ylab="Proportion")
dev.off()

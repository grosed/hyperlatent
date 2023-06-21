library('plyr')
library('HyperG') 
library('igraph')

library(praxi)

#source('subsample.R')
#source('../../mcmc/graphsummary.R') 

######################################################################
# get DBLP data from  https://xueyumao.github.io/coauthorship.html : 
# download https://xueyumao.github.io/data/coauthorship/DBLP_bipartite.zip
# into current folder
######################################################################

## load the community memberships 
dblp_comms = read.table(file="./DBLP_bipartite/DBLP3_bipartite_community.txt")
sum(dblp_comms$V2 %in% c(4,5,6))  
sum(dblp_comms$V2 %in% c(1,2,3)) 

## load the incidence matrix (note: there are 5 of these, numbered 1-5)
dblp = read.table(file="./DBLP_bipartite/DBLP3_bipartite_adjacency.txt")
nc1 = length(unique(dblp$V1)) 
nc2 = length(unique(dblp$V2)) 
colnames(dblp) = c("paperID", "authorID") 

## create an edgelist
auths = unique(dblp$authorID)
N = length(auths)
papers = unique(dblp$paperID)
elst = vector(mode="list", length=length(papers))
for (e in 1:length(elst)){
    elst[[e]] = dblp$authorID[ dblp$paperID == papers[e] ]
}

## check density/other summaries
N = length(unique(unlist(elst))) 
elen = unlist(lapply(elst, length))
K = max(elen) # 24 authors on a paper
hist(elen) # looks quite reasonable

## subset so that hyperedges are 'right' size
Kmax = 5
e_kp = (elen>1) & (elen<Kmax)
elst = elst[e_kp]
elst = lapply(elst, sort)
elst = unique(elst)

## create hypergraph object (useful for subsetting)
dblp_h = hypergraph_from_edgelist(elst) # take a little time

## check if connected
hypergraph.is.connected(dblp_h) # this gives a false

## isolate largest connected component - contains 9798 nodes
conn = hyp_conn_comp(dblp_h)
sub = which(conn$membership == which(conn$csize==max(conn$csize)))
dblp_h_lcc = induced_hypergraph(dblp_h, sub) 
dblp_elst_lcc = hypergraph_as_edgelist(dblp_h_lcc) # NOTE: nodes are labelled as CHARACTERS

## sort elst and make sure no repeated hyperedges
dblp_elst_lcc = lapply(dblp_elst_lcc, as.numeric) # convert to numbers
dblp_elst_lcc = lapply(dblp_elst_lcc, sort)
dblp_elst_lcc = unique(dblp_elst_lcc)
dblp_elst_lcc = lapply(dblp_elst_lcc, as.character) # convert back to characters

######################################################################
######## random walk subsample to assess 'additional nodes' ##########
######################################################################

## take larger subsample, and 'hold out' nodes on the periphery
degs = sort(hdegree(dblp_h_lcc))

set.seed(843)
nd_kp = random_walk_subsample(dblp_elst_lcc, 38, c_ret = 0.15, nmin = 90)  
h_dblp_rw = induced_hypergraph(dblp_h_lcc, nd_kp) # 41 nodes
elst_dblp_rw = hypergraph_as_edgelist(h_dblp_rw)

plot(h_dblp_rw) 

## make sure edges are unique (some may be repeated)
elst_dblp_rw = lapply(elst_dblp_rw, as.numeric) # convert to numbers
elst_dblp_rw = lapply(elst_dblp_rw, sort)
elst_dblp_rw = unique(elst_dblp_rw)
elst_dblp_rw = lapply(elst_dblp_rw, as.character) # convert back to characters

## identify nodes on the periphery
g_dblp_sb = as.graph(h_dblp_rw)
cent = eigen_centrality(g_dblp_sb)
wv_add = names(cent$vec)[match( sort(cent$vec)[1:9], cent$vec)] # these are the 'new' indices

##### relabel everything #####

## separate into subset to fit on and subset to calc predictives for
wv_kp = setdiff(unique(unlist(elst_dblp_rw)), wv_add)
wv_full = c(wv_add, wv_kp) # 106 in total

## 'full' subsample
dblp_rw_full = induced_hypergraph(dblp_h_lcc, wv_full) 
rw_full_elst = lapply(hypergraph_as_edgelist(dblp_rw_full), as.numeric)
rw_full_elst = lapply(rw_full_elst, sort) # sort hyperdges
rw_full_elst = unique(rw_full_elst) # remove repeated hyperedges

## take subsample to fit on by finding hyperedges which only contain indices wv_kp
id_rw_sub = unlist(lapply( rw_full_elst, function(x){ all( x %in% wv_kp) } ))
rw_sub_elst = rw_full_elst[id_rw_sub]
rw_sub_elst = lapply(rw_sub_elst, sort) # sort hyperdges
rw_sub_elst = unique(rw_sub_elst) # remove repeated hyperedges

## make a reference list for all vertices to relabel as 1:N(str)
## NOTE: not the same as wv_kp and wv_full since some 'new' nodes only interact with 'old' nodes
srt_kp = sort(unique(unlist(rw_sub_elst))) 
srt_wv = sort(setdiff(unique(unlist(rw_full_elst)), unique(unlist(rw_sub_elst)))) 
srt_full = as.numeric(c(srt_kp, srt_wv))

## relabel indices to be 1:N(str)
rw_sub_elst = lapply(rw_sub_elst, function(x){match(x, srt_full)} )
rw_sub_elst = lapply(rw_sub_elst, sort) # ensure e's i<j<...
rw_full_elst = lapply(rw_full_elst, function(x){match(x, srt_full)} )
rw_full_elst = lapply(rw_full_elst, sort) # ensure e's i<j<...

## plot for sanity
par(mfrow=c(1,2))
plot(hypergraph_from_edgelist(rw_sub_elst))
plot(hypergraph_from_edgelist(rw_full_elst))

N = length(unique(unlist(rw_sub_elst))) # may differ from length of wv_kp
Nstr = length(unique(unlist(rw_full_elst))) - N # may differ from length of wv_add

## save output
id_new = unlist(lapply(rw_full_elst, function(x){sum(unlist(lapply(rw_sub_elst, function(y){ identical(x, y) } ) ) ) } )) 
rw_extra_elst = rw_full_elst[id_new==0]

save(rw_sub_elst, file="dblp_rw_to_fit.RData") # subset to fit on
save(rw_full_elst, file="dblp_rw_full.RData") # full data
save(rw_extra_elst, file="dblp_rw_extra.RData") # 'added' hyperedges

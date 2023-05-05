library('HyperG') 
library('igraph')
source('subsample.R')

######################################################################
############## open grocery data from arules #########################
######################################################################

library(arules)
data(Groceries)
gro = hypergraph_from_incidence_matrix( t(Groceries@data) ) # is the data matrix
dim(gro$M) # 9835 (transactions) x 169 (categories)

## take nodes to be categories (ie items) and hyperedges to be transactions
nd_deg = colSums(gro$M) # tells me the degree of each node
e_deg = rowSums(gro$M) # tells me the size of each hyperedge

## remove hyperedges with deg < 2, > Kmax
Kmax = 8
mat = gro$M
e_deg = rowSums(mat)
e_kp = (e_deg > 1) & (e_deg < Kmax)
mat = mat[e_kp,]

## create an edgelist
elst = lapply( 1:dim(mat)[1], function(x){ as.numeric(names(which(mat[x,]))) } )
elst = lapply( elst, sort )
elst = unique( elst ) # remove multiple copies of hyperedges
h_gro = hypergraph_from_edgelist(elst)
elst_gro = hypergraph_as_edgelist(h_gro) # store hyperedges as characters

## sanity check: one connected component - if not, need to isolate largest connected component
elen = unlist( lapply( elst, length ) )
cc = hyp_conn_comp( h_gro ) # gives a single connected component

######################################################################
###################### random subsample ##############################
######################################################################

N = length(unique(unlist(elst))) # 167 nodes in total

set.seed(98745)
nsml = 46
n_kp = sample(1:N, nsml)
h_rndm = induced_hypergraph(h_gro, n_kp)
elst_rndm = hypergraph_as_edgelist(h_rndm)
K_rndm = max(unlist(lapply(elst_rndm, length)))

## make sure edges are unique (some may be repeated)
elst_rndm = lapply(elst_rndm, as.numeric) # convert to numbers
elst_rndm = lapply(elst_rndm, sort)
elst_rndm = unique(elst_rndm)
elst_rndm = lapply(elst_rndm, as.character) # convert back to characters

## identify nodes on the periphery
g_gro_sb = as.graph(h_rndm)
cent = eigen_centrality(g_gro_sb)
wv_add = names(cent$vec)[match( sort(cent$vec)[1:5], cent$vec)] # these are the 'new' indices

##### relabel everything #####

## separate into subset to fit on and subset to calc predictives for
wv_kp = setdiff(unique(unlist(elst_rndm)), wv_add)
wv_full = c(wv_add, wv_kp) # 106 in total

## 'full' subsample
gro_rndm_full = induced_hypergraph(h_gro, wv_full) 
rndm_full_elst = lapply(hypergraph_as_edgelist(gro_rndm_full), as.numeric)
rndm_full_elst = lapply(rndm_full_elst, sort) # sort hyperdges
rndm_full_elst = unique(rndm_full_elst) # remove repeated hyperedges

## take subsample to fit on by finding hyperedges which only contain indices wv_kp
id_rndm_sub = unlist(lapply( rndm_full_elst, function(x){ all( x %in% wv_kp) } ))
rndm_sub_elst = rndm_full_elst[id_rndm_sub]
rndm_sub_elst = lapply(rndm_sub_elst, sort) # sort hyperdges
rndm_sub_elst = unique(rndm_sub_elst) # remove repeated hyperedges

## make a reference list for all vertices to relabel as 1:N(str)
srt_kp = sort(unique(unlist(rndm_sub_elst))) 
srt_wv = sort(setdiff(unique(unlist(rndm_full_elst)), unique(unlist(rndm_sub_elst)))) 
srt_full = as.numeric(c(srt_kp, srt_wv))

## relabel indices to be 1:N(str)
rndm_sub_elst = lapply(rndm_sub_elst, function(x){match(x, srt_full)} )
rndm_sub_elst = lapply(rndm_sub_elst, sort) # ensure e's i<j<...
rndm_full_elst = lapply(rndm_full_elst, function(x){match(x, srt_full)} )
rndm_full_elst = lapply(rndm_full_elst, sort) # ensure e's i<j<...

## plot for sanity
par(mfrow=c(1,2))
plot(hypergraph_from_edgelist(rndm_sub_elst))
plot(hypergraph_from_edgelist(rndm_full_elst))

N = length(unique(unlist(rndm_sub_elst))) # may differ from length of wv_kp
Nstr = length(unique(unlist(rndm_full_elst))) - N # may differ from length of wv_add

## save output
id_new = unlist(lapply(rndm_full_elst, function(x){sum(unlist(lapply(rndm_sub_elst, function(y){ identical(x, y) } ) ) ) } )) 
rndm_extra_elst = rndm_full_elst[id_new==0]

save(rndm_sub_elst, file="gro_rndm_to_fit.RData") # subset to fit on
save(rndm_full_elst, file="gro_rndm_full.RData") # full data
save(rndm_extra_elst, file="gro_rndm_extra.RData") # 'added' hyperedges

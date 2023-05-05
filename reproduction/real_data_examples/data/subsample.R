library(snowboot)
library(igraph)
library(HyperG)

hypelst_to_graph <- function(elst, wght=TRUE){
    nd_ids = sort(unique(unlist(elst)))
    N = length(nd_ids)
    adj = matrix(0, N, N)
    for (i in 1:(N-1)){
        for (j in (i+1):N){
            e = sum(unlist(lapply(elst, function(x){ sum(c(nd_ids[i], nd_ids[j]) %in% x) == 2})))
            if (wght==TRUE){
                adj[i,j] = e
                adj[j,i] = e
            } else if (wght==FALSE){
                adj[i,j] = e > 0
                adj[j,i] = e > 0
            }
        }
    }
    return( adj )
}

hyp_conn_comp <- function(h){ # find connected components
    g = graph_from_adjacency_matrix(as.matrix(hypergraph_as_adjacency_matrix(h)), mode="undirected")
    comps = components(g)
    return( comps )
}

hyp_kcore <- function(h){ # find k core - VERY slow
    g = graph_from_adjacency_matrix(as.matrix(hypergraph_as_adjacency_matrix(h)), mode="undirected")
    core = coreness(g)
    kvals = unique(core)
    return( comps )
}

snowball_hypergraph_indices <- function(elst, input_seed = NA, nw = 2){
    ## elst is an edge list for a hypergraph

    ## create hypergraph
    h = hypergraph_from_edgelist(elst)
    ## convert to graph and define network object
    g = as.graph(h)
    degs = degree(g)
    N = vcount(g)
    net = list( edges = as_edgelist(g), degree = degs, n = N )

    ## snowball sample the graph
    if( is.na(input_seed) == TRUE ){
        seed = sample(1:N, 1, prob=degs/sum(degs))
    } else {
        seed = input_seed
    }
    sb = sample_about_one_seed(net, seed=seed, n.wave=nw )

    ## get list of newly added vertices as each wave
    sb_na = lapply(sb, unique)
    for (i in 2:length(sb_na)){
        sb_na[[i]] = setdiff( sb_na[[i]], unique(unlist(sb_na[1:(i-1)])) )
    }

    return( list(sb=sb, sb_na=sb_na, seed=seed) )
}

sample_hyperedge_node <- function(elst, nd){
    ## extract hyperedges incident to nd
    es_adj_nd = elst[unlist(lapply(elst, function(x){ nd %in% x } ))]
    ## sample uniformly from hyperedges
    n_es_adj = length(es_adj_nd)
    smp_e = unlist(es_adj_nd[sample(1:n_es_adj, 1)])
    return(smp_e)
}

random_walk_subsample <- function(elst, seed, c_ret = 0.15, nmin = 50){
                                        # elst is an edge list for a hypergraph
                                        # seed is starting node
                                        # wght = T/F, whether projected graph is weights
                                        # c_ret is probability of returning to seed node
                                        # nmin is minimum length of newly sampled nodes

   
    ## determine random walk subsample
    nsmp = seed # start with seed node, append additional nodes later on

    while( length(nsmp) <= nmin ){
        ## sample random walk, at each stage can go back to the start with probability c_ret

        ## start at seed node, and set 'return index' to 0
        cur_seed = seed
        ret = 0 
        while( ret == 0 ){
            new_nds = sample_hyperedge_node(elst, cur_seed) 
            new_nds = new_nds[ new_nds != cur_seed ] # remove current seed node (included by default)
            ret = rbinom(1, 1, c_ret)
            cur_seed = new_nds[sample(length(new_nds), 1)]
            nsmp = append(nsmp, new_nds)
        }
        nsmp = unique(nsmp)
        print( length(nsmp) )
    }

    return( nsmp )
}

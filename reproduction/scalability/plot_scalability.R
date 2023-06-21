## plot output from scalability study
load(file="./output/scalability_dens.RData") # dens: smaller/larger, reps, density, nodes
load(file="./output/scalability_nmat.RData") # nmat: rsmaller/larger, reps, nvalues
load(file="./output/scalability_times.RData") # timings: r smaller/large, reps, DA/non, nodes

Nvec = c(10, 30, 50, 100) # approx N values
nRep = dim(nmat)[2]

## plot densities
pdf(file="./plots/scalability_sim_density.pdf", width=7, height=7)
par(mfrow=c(2,2))
## order 2 densities for each N, r_sml
plot( nmat[1,1,], dens[1,1,1,] , ylim=c(0, max(dens[,,1,], na.rm=TRUE) ), main="(R1)", xlab="N", ylab="Density of k=2 hyperedges", pch=19 )
for( rep in 2:nRep ){ points( nmat[1,rep,], dens[1,rep,1,] , pch=19 ) }
## order 2 densities for each N, r_lrg
plot( nmat[2,1,], dens[2,1,1,] , ylim=c(0, max(dens[,,1,], na.rm=TRUE) ), main="(R2)", xlab="N", ylab="Density of k=2 hyperedges" , pch=19 )
for( rep in 2:nRep ){ points( nmat[2,rep,], dens[2,rep,1,] , pch=19 ) }

## order 3 densities for each N, r_sml
plot( nmat[1,1,], dens[1,1,2,] , ylim=c(0, max(dens[,,2,], na.rm=TRUE) ), main="(R1)", xlab="N", ylab="Density of k=3 hyperedges" , pch=19 )
for( rep in 2:nRep ){ points( nmat[1,rep,], dens[1,rep,2,] , pch=19 ) }
## order 3 densities for each N, r_lrg
plot( nmat[2,1,], dens[2,1,2,] , ylim=c(0, max(dens[,,2,], na.rm=TRUE) ), main="(R2)", xlab="N", ylab="Density of k=3 hyperedges" , pch=19 )
for( rep in 2:nRep ){ points( nmat[2,rep,], dens[2,rep,2,] , pch=19 ) }
dev.off()

## plot relative timings
pdf(file="./plots/scalability_sim_R1_and_R2.pdf", width=8, height=4)
par(mfrow=c(1,2))
## raw timings for DA scheme
plot( nmat[1,1,], timings[1,1,1,], ylim = c(0, max(timings[,,1,], na.rm=TRUE)), ylab="DA average per iteration cost", xlab="N", pch=19, xlim=c(0,105))
for ( rep in 2:nRep){ points( nmat[1,rep,], timings[1,rep,1,], pch=19 ) }
for ( rep in 1:nRep){ points( nmat[2,rep,], timings[2,rep,1,], pch=17, col='red' ) }
##ratio of DA timings/non DA timings
ymx = max(timings[,,1,]/timings[,,2,], na.rm=TRUE)
plot( nmat[1,1,], timings[1,1,1,]/timings[1,1,2,], ylim = c(0, ymx), ylab="DA / non-DA", xlab="N", pch=19, xlim=c(0,105))
for ( rep in 2:nRep){ points( nmat[1,rep,], timings[1,rep,1,]/timings[1,rep,2,], pch=19 ) }
for ( rep in 1:nRep){ points( nmat[2,rep,], timings[2,rep,1,]/timings[2,rep,2,], pch=17, col='red' ) }
legend("topright", legend=c("(R1)", "(R2)"), col=c('black', 'red'), pch = c(19,17) )
dev.off()

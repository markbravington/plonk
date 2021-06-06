## Pinocchio's Dolphin (Delfinus mendax) from Bilateral Bay: L- and R-handed

# Lethal sampling; adult age UNKNOWN, juve age OK; "mammals"...
# ... ie knife-edge maturity; all adults equal WRTO reprod & survival; constant adult survival
# "Sampling season" is after "breeding season"
# Only female samples used; only female pop dyn considered; MOPs are actually MDPs
# MOPs only (HSPs not used yet)

library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun

### ACTION STARTS HERE

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))
print( load( 'samp_delfi_D_raw.rda')) # object samp_notog1
print( names( samp_delfi_D))
print( head( samp_delfi_D$Samps))

samp_delfi_D <- "READ THE DAMN SCRIPT AND EDIT APPROPRIATELY; DON'T JUST RUN IT!"

# Summarize sampling & kin & other data for delfi_D, and teach the lglk function about it
lglk_delfi_D <- generic_lglk_MOP_fuzzall_mammal # a copy that will be taught about delfi_D data
env <- boring_data_prep_delfi_D( samp_delfi_D, prev_env=environment( lglk_delfi_D))
environment( lglk_delfi_D) <- env # now, lglk_<blah> will always know eg what...
# ... catches, MOPs, comparisons etc were, without having to pass that stuff in every time it's called

# Fitting is sloooooow, so...
cat( 'Fitting is sloooow, so gonna use previously-fitted params... *hopefully* for this exact dataset\n')
lglk_delfi_D( c( 9.11279 0.02291) )
if( FALSE) {
  starto <- c( log( 100000), 0.01)
  lglk_delfi_D( starto) # just check it doesn't crash! NB slow.
  fitto <- nlminb( starto, NEG( lglk_delfi_D), control=list( trace=1))
}

# Howzit?
Amat <- env$Amat
Brange <- with( env, y0 %upto% ylast) # what years were estimated?

tru_Nadf <- sumover( samp_delfi_D@secret$N_sya[SLICE='F', Brange, Amat:40], 'A')

maxyplot <- max( c( tru_Nadf, env$N)) * 1.05
plot( Brange, tru_Nadf, type='l', ylim=c( 0, maxyplot), ylab='Nadf')
points( Brange, env$N, col='green', pch=16)

#######
# Insert your own diagnostic code here--- eg to look at ppns by LL LR etc
# Look at fit_delfi_A.r
# first compute Expected[MOPs]
# then sum it up in an informative way
# and do the same for observed MOPs (n_MOP[...])
# and compare them

######################
# A WRONG WAY... one of many, no doubt
# Actually even this is harder than i thought so I gave up
# Basically you get "impossible" POPs ( 

# Bozo et al (2020b), just using measured age ("it's unbiased after all")
# Widely cited
# ie "Equations are for ivory-tower pointy-heads The Reality Is I'm a Practical Man"
# But see Shaw, GB on the subject of "Practical Men" (I think GBS, but Google is unforthcoming)

# Use ideal mammal (Acme Archipelago code) but keeping Amat
lglk_delfi_D_ignore <- generic_lglk_MOP_ideal_mammal
samp_delfi_D_ignore <- samp_delfi_D
samp_delfi_D_ignore$Samps <- within( samp_delfi_D_ignore$Samps, {
  A <- Fuzzage
  Fuzzage <- NULL
  poss_par <- A >= Amat
  poss_off <- A < Amat
})
swap_POP <- with( samp_delfi_D_ignore, with( Samps, {
  poss_off[ POPs[,1]] & poss_par[ POPs[,2]] # "wrong way round"
}))
samp_delfi_D_ignore$POPs[ swap_POP,] <- samp_delfi_D_ignore$POPs[ swap_POP,2:1]

drop_POP <- with( samp_delfi_D_ignore, with( Samps, {
  !poss_par[ POPs[,1]] | !poss_off[ POPs[,2]]
}))
samp_delfi_D_ignore$POPs <- samp_delfi_D_ignore$POPs[ !drop_POP,]

# Also a few "impossibles" on maturity grounds
drop_POP <- with( samp_delfi_D_ignore, with( Samps, {
  Yad <- Y[ POPs[,1]]
  Aad <- A[ POPs[,1]]
  Bju <- Y[ POPs[,2]] - A[ POPs[,2]]
  Bad <- Yad - Aad
  # Code from gen_lglk_IDEAL:
  believable_at_face_value  <-     
	    ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
	    ( Bju >= Amat + Bad )   # was ad mature at B[ju] ?
	!believable_at_face_value
}))
samp_delfi_D_ignore$POPs <- samp_delfi_D_ignore$POPs[ !drop_POP,]

# Now treat whatever's left like exact data
env2 <- boring_data_prep_delfi_A( samp_delfi_D_ignore, prev_env=environment( lglk_delfi_D_ignore))
environment( lglk_delfi_D_ignore) <- env2

starto <- c( log( 100000), 0.01)

lglk_delfi_D_ignore( starto) # just check!
fitto2 <- nlminb( starto, NEG( lglk_delfi_D_ignore))

maxyplot <- max( c( maxyplot, env2$N)) * 1.05
plot( Brange, tru_Nadf, type='l', ylim=c( 0, maxyplot), ylab='Nadf')
points( Brange, env$N, col='green', pch=16)

# Brange may not be the same after Bozo et al have finished with it
with( env2, points(dimseq(N), N, col = "orange", pch = 16))
# NB dimseq returns a vector with 1D arrays, rather than a length-1 list...
# ... unless drop=F

## TBF there are lots of other WRONG ways to analyse this Fuzzage dataset--- 
## feeeeeeeeel freeeeeee to try them...

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
print( load( 'samp_delfi_B_raw.rda')) # object samp_notog1
print( names( samp_delfi_B))
print( head( samp_delfi_B$Samps))

samp_delfi_B <- "READ THE DAMN SCRIPT AND EDIT APPROPRIATELY; DON'T JUST RUN IT"

# Summarize sampling & kin & other data for delfi_B, and teach the lglk function about it
lglk_delfi_B <- generic_lglk_MOP_LR_mammal # a copy that will be taught about delfi_B data
env <- boring_data_prep_delfi_B( samp_delfi_B, prev_env=environment( lglk_delfi_B)) 
environment( lglk_delfi_B) <- env # now, lglk_<blah> will always know eg what...
# ... catches, MOPs, comparisons etc were, without having to pass that stuff in every time it's called

starto <- c( log( 100000), 0.01, logit( 0.4))
lglk_delfi_B( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_delfi_B))

# Howzit?
Amat <- env$Amat
Brange <- dimseq( env$N) # what years get estimated?

tru_Nadf <- sumover( samp_delfi_B@secret$N_sya[SLICE='F', Brange, Amat:40], 'A')

maxyplot <- max( c( tru_Nadf, env$N))
plot( Brange, tru_Nadf, type='l', ylim=c( 0, maxyplot))
points( Brange, env$N, col='green', pch=16)

#######
# Insert your own diagnostic code here--- eg to look at ppns by LL LR etc
# Look at fit_delfi_A.r
# first compute Expected[MOPs]
# then sum it up in an informative way
# and do the same for observed MOPs (n_MOP[...])
# and compare them

######################
# THE WRONG WAY
# Following Bozo et al (2020), applying "noHands" model that ignores Handedness 
# ie "I just want a quick answer equations are for ivory-tower pointy-heads The Reality Is I'm a Practical Man"
# See Shaw, GB on the subject of "Practical Men"
lglk_delfi_B_noHands <- generic_lglk_MOP_ideal_mammal
env2 <- boring_data_prep_delfi_B( samp_delfi_B, prev_env=environment( lglk_delfi_B_noHands)) 
# ... fudge to use the noHands data summaries. NB this works directly on environment( lglk...) !
env2$n_MOP <- env2$n_MOP_noHands
env2$n_comp_MOP <- env2$n_comp_MOP_noHands
environment( lglk_delfi_B_noHands) <- env2

starto <- c( log( 100000), 0.01)
lglk_delfi_B_noHands( starto) # just check!
fitto2 <- nlminb( starto, NEG( lglk_delfi_B_noHands))

points(Brange, env2$N, col = "orange", pch = 16)
# Fairly bad bias... though I really had to crank up the effect!

##############
# Diagnostics...
# if you actually HAVE the L/R and just didn't use them, then 
# well-chosen Obs&Exp diagnostics should reveal a problem


##########
# Interesting stuff that I didn't look at:
# what if you didn't trust the Dominant/Recessive model?
# but what if you did have external estimate of Pr[R/L], eg from strandings (presumably unbiased WRTO L/R)?


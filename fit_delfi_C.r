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
print( load( 'samp_delfi_C_raw.rda')) # object samp_notog1
print( names( samp_delfi_C))
print( head( samp_delfi_C$Samps))

samp_delfi_C <- "READ THE DAMN SCRIPT AND EDIT APPROPRIATELY; DON'T JUST RUN IT!"

# Summarize sampling & kin & other data for delfi_C, and teach the lglk function about it
lglk_delfi_C <- generic_lglk_MOP_fuzzage_adult_mammal # a copy that will be taught about delfi_C data
env <- boring_data_prep_delfi_C( samp_delfi_C, prev_env=environment( lglk_delfi_C))
environment( lglk_delfi_C) <- env # now, lglk_<blah> will always know eg what...
# ... catches, MOPs, comparisons etc were, without having to pass that stuff in every time it's called

# Have a look at the deconvolution stuff:
matplot( unclass( env$samp_Fuzzad), type='l', main='Observed *adult* age from P.E.', xlab='Age', ylab='Freq')
matplot( unclass( env$est_Pr_A_SampY), type='l', main='Est ad age after deconvo--- add 5 to this!', xlab='Age', ylab='Freq')

starto <- c( log( 100000), 0.01)
lglk_delfi_C( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_delfi_C))

# Howzit?
Amat <- env$Amat
Brange <- dimseq( env$N) # what years get estimated?

tru_Nadf <- sumover( samp_delfi_C@secret$N_sya[SLICE='F', Brange, Amat:40], 'A')

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
# A WRONG WAY
# Bozo et al (2020b), just using measured age ("it's unbiased after all")
# Widely cited
# ie "Equations are for ivory-tower pointy-heads The Reality Is I'm a Practical Man"
# But see Shaw, GB on the subject of "Practical Men"


# Use ideal mammal (Acme Archipelago code) but keeping Amat
lglk_delfi_C_ignore <- generic_lglk_MOP_ideal_mammal
env2 <- boring_data_prep_delfi_C( samp_delfi_C, prev_env=environment( lglk_delfi_C_ignore))
# ... fudge to change the names
evalq( envir=env2, { # will happen "live" inside env2
  names( dimnames( n_MOP))[3] <- 'A' # instead of "Fuzzage"
  names( dimnames( n_comp_MOP))[3] <- 'A'
  Aad_range <- Fuzzad_range
  
  # Better zero MOPs & comps for "impossibles", following Bozo et al (2020b)
  OK <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=Aad_range, {
    Bad <- Yad - Aad
    
    # 1/0 value:
    1 * ( Bju <= Yad) *          # was ad still alive at B[ju] ?
        ( Bju >= Amat + Bad)    # was ad mature at B[ju] ?
  })
  
  n_MOP <- n_MOP * OK
  n_comp_MOP <- n_comp_MOP * OK
})
 
 
environment( lglk_delfi_C_ignore) <- env2

starto <- c( log( 100000), 0.01)
lglk_delfi_C_ignore( starto) # just check!
fitto2 <- nlminb( starto, NEG( lglk_delfi_C_ignore))

maxyplot <- max( c( maxyplot, env2$N)) * 1.05
plot( Brange, tru_Nadf, type='l', ylim=c( 0, maxyplot), ylab='Nadf')
points( Brange, env$N, col='green', pch=16)
points(Brange, env2$N, col = "orange", pch = 16)

## TBF there are lots of other WRONG ways to analyse this Fuzzage dataset--- 
## feel freeeeeee to try them...

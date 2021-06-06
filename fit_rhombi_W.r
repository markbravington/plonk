## Cobalt Squarehad (Rhombichthys cyanorosea) from Whyalla

library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun

### ACTION STARTS HERE

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))
print( load( 'samp_rhombi_W_raw.rda')) # object samp_notog1
print( names( samp_rhombi_W))
print( head( samp_rhombi_W$Samps))
scatn( 'Rhombi W #POPs = %i', nrow( samp_rhombi_W$POPs))

samp_rhombi_W <- "READ THE DAMN SCRIPT AND EDIT APPROPRIATELY; DON'T JUST RUN IT"

# Summarize sampling & kin & other data for rhombi_W, and teach the lglk function about it
lglk_rhombi_W <- generic_lglk_Covar # a copy that will be taught about rhombi_W data
env <- boring_data_prep_rhombi_W( samp_rhombi_W, prev_env=environment( lglk_rhombi_W))
environment( lglk_rhombi_W) <- env # now, lglk_<blah> will always know about data

starto <- c( rep( log( 100000), 2), 0.5) # 
lglk_rhombi_W( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_rhombi_W))

# How did  we go?
cat( 'Abund:\n')
rbind( Est=env$Nad, Tru=samp_rhombi_W@secret$Nad)

cat( 'Fec\n')
rbind( Est=with( env, fec_Co_m / fec_Co_m[3]), 
    Tru=with( samp_rhombi_W@secret, fec_Co / fec_Co[3])
  )

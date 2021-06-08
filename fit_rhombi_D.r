## Cobalt Squarehad (Rhombichthys cyanorosea) from Derwent

library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun

### ACTION STARTS HERE

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))
print( load( 'samp_rhombi_D_raw.rda')) # object samp_notog1
print( names( samp_rhombi_D))
print( head( samp_rhombi_D$Samps))
scatn( 'Rhombi D #POPs = %i', nrow( samp_rhombi_D$POPs))

samp_rhombi_D <- "READ THE DAMN SCRIPT AND EDIT APPROPRIATELY; DON'T JUST RUN IT"

# Summarize sampling & kin & other data for rhombi_D, and teach the lglk function about it
lglk_rhombi_D <- generic_lglk_cartoon # a copy that will be taught about rhombi_D data
env <- boring_data_prep_rhombi_D( samp_rhombi_D, prev_env=environment( lglk_rhombi_D))
environment( lglk_rhombi_D) <- env # now, lglk_<blah> will always know about data

starto <- rep( log( 100000), 2) # 
lglk_rhombi_D( starto) # just check it doesn't crash!

fitto <- nlminb( starto, NEG( lglk_rhombi_D))

# How did  we go?
cat( 'Est:\n')
env$Nad # estimated
cat( 'Tru:\n')
samp_rhombi_D@secret$Nad # trooo

# Compare with bleedin' obvious...
Nobvious <- with( env, n_comp_POP / n_POP)

# How would aggregate version go?
Nsexlumpo <- with( env, 2 * sum( n_comp_POP) / sum( n_POP))
sum( Nobvious)
# ... pretty close, as you'd hope with nonselective sampling

# For interest's sake, is there any link to Co?
# Females?
with( samp_rhombi_D, with( Samps, {
  table( Co[ Sex=='F'])
}))
# All female adults have Co==1 , so nothing to see there !

Co_m_samp <- with( samp_rhombi_D, with( Samps, {
  table( Co[ Sex=='M'])
}))

Co_m_POP <- with( samp_rhombi_D, with( Samps, {
  table( Co[ POPs[,1][ Sex[ POPs[,1]]=='M'] ])
}))

Co_m_POP / Co_m_samp
# Hmm! Doesn't affect the answer, though. And that is because <complete in your own words>

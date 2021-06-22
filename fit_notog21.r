## Simplest conceivable "general" CKMR model:

# Lethal sampling; adult age UNKNOWN, juve age OK; "mammals"...
# ... ie knife-edge maturity; all adults equal WRTO reprod & survival; constant adult survival
# "Sampling season" is after "breeding season"
# Only female samples used; only female pop dyn considered; POPs are actually MDPs
# POPs only (HSPs not used yet)

## NB This is 2021 update of 2019 Halifax code--- may not follow _exactly_ the same
## pattern as the other 2021 examples
## Caveat lector!

library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( debug)     # so we can watch the fun

# source( 'borings.r')
# source( 'lglks.r')

### ACTION STARTS HERE

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))
print( load( 'samp_notog1.rda')) # object samp_notog1
print( names( samp_notog1))
print( head( samp_notog1$Samps))

lglk_notog <- generic_lglk_fish # a copy that will be taught about delfi_C data
env <- boring_data_prep_fish( samp_notog1, prev_env=environment( lglk_notog))
environment( lglk_notog) <- env # now, lglk_notog will always know eg what...
# ... catches were, etc, without having to pass that stuff in every time it's called

## Where shall we start..?
# Cheat this time, based on true numbers in population...
# starting values don't affect final estimate (in principle) so...
# this is not an important cheat !
param_start <- cheat_ruff_tru_pars( samp_notog1, env)

# env$verbose <- TRUE # debugging only
print( lglk_notog( param_start))

# Used optim() rather than nlminb() here; I'm never sure which to use
# (for real data, I often try both) but with optim() you can just say hessian=TRUE et voila!
# NB do always use method='BFGS' or 'L-BFGS-B' with optim()
fit <- optim(
    param_start,
    NEG( lglk_notog),
    method='BFGS',
    hessian=TRUE
    # ,control=list( trace=5)
  )

print( fit) # things OK?
print( exp( fit$par)) # on unlogged scale
lglk_notog( fit$par) # one final call, to ensure things in env are "best estimates"

## What was the estimated pop dyn?
cat( 'EST numbers-at-year-and-age (FEMALES):\n')
print( env$N_ya)

cat( 'TRUE ditto:\n')
print( attr( samp_notog1, "secret")$N_sya[ SLICE='F', 28:40, 4:12])

# ... and you will find other interesting estimates in env, too

# How many kin-pairs did we expect? I'm not sure if this should match exactly for both POPs and HSPs with fish,
# because the pop dyn model is not exactly correct (eg exactly geometric numbers-at-age in year "0")...
# ... but it seems to!

with( env,
   scatn( 'Obs & Exp POPs: %5.0f, %5.2f', sum( n_POP),
       sum( Pr_POP * n_comp_POP)))
with( env,
   scatn( 'Obs & Exp HSPs: %5.0f, %5.2f', sum( n_HSP),
       sum( Pr_HSP * n_comp_HSP)))


# More fancy diagnostics: compare E[POP] & Obs[POP] by some split of categories

# Some things to look at:
# ERRO on different scales eg equiv number of 5yo, or of 10yo

scatn( "Now look at 'notog_summaries.r' for SPRR etc")



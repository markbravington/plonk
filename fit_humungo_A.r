# Lethal sampling; HSP s
# ... ie knife-edge maturity; all adults equal WRTO reprod & survival; constant adult survival
# "Sampling season" is after "breeding season"
# Only female samples used; only female pop dyn considered; MHSPs only (i hope)


### ACTION STARTS HERE
library( mvbutils)  # various low-level utilities, eg cq()
library( offarray)  # arrays/vectors don't have to start at 1 any more! Freedom!
library( atease)    # x@a <=> attr( x, "a")
library( debug)     # so we can watch the fun

## Load and set up the data
od <- options( digits=7)
on.exit( options( od))

print( load( 'samp_humungo_A_raw.rda'))
print( names( samp_humungo_A))
print( head( samp_humungo_A$Samps))
print( head( samp_humungo_A$MHSPs))
print( ls( samp_humungo_A@public))

# A, Y: Age-at-sampling, Year-of-sampling
# poss_off: ?use as candidate HSP (or O in POP) ?
# poss_par: ?use as candidate P in POP?
# Amat: (female) maturity --- in @public


# eg comparison-order of pairs, etc, yawn zzzz

lglk_humungo_A <- generic_lglk_HSP_mammal # a copy that will be taught about humungo_A data
env <- boring_data_prep_humungo_A( samp_humungo_A, prev_env=environment( lglk_humungo_A))
environment( lglk_humungo_A) <- env # now, lglk_<blah> will always know eg what...
# ... catches, MOPs, comparisons etc were, without having to pass that stuff in every time it's called


## Where shall we start..?
param_start <- c( log_N0=log( 1000), RoI=0.01, log_Z=log( 0.2))

print( lglk_humungo_A( param_start))


# Use optim() instead for nlminb(), just for a change
fit_hsp <- optim(
    param_start,
    NEG( lglk_humungo_A),
    method='BFGS',
    hessian=TRUE
  )
print( fit_hsp) # things OK?

lglk_humungo_A( fit_hsp$par) # one more run at best pars, to ensure env has The Right Stuff

## What was the estimated pop dyn?

scatn( 'RoI and Z estimate: %5.2f %5.2f', env$RoI, env$Z)

cat( 'Nfad est and truth:\n')
Aad_range <- env$Amat %upto% dimrange( samp_humungo_A@secret$N_sya)['A',2]
trooth <- with( env, rowSums( samp_humungo_A@secret$N_sya[ SLICE='F', Bju_range, Aad_range]))  # num Fem ads
print( rbind( EST_Nfad=with( env, N[ Bju_range]), TRU_Nfad= trooth))

# matplot( t( nest_tru), type='l', ylim=c( 0, max( c( nest_tru))))

# How many HSPs did we expect? This should match exactly since "N" is a free param, so it's just a sanity check
with( env,
   scatn( 'Obs & Exp MHSPs: %5.0f, %5.2f', round( sum( n_MHSP)), sum( Pr_MHSP * n_comp_MHSP)))


# To demonstrate variance calculation, let's look at *final* abundance
log_final_abund <- function( params) {
    lglk_humungo_A( params) # set up pop dyn
    envo <- environment( lglk_humungo_A) # paranoid programming--- env
    Ylast <- lastel( envo$N)
    log_N_final <- log( envo$N[ Ylast])
    names( log_N_final) <- sprintf( 'Y%i', Ylast)
  return( log_N_final)
  }

# check...
exp( log_final_abund( fit_hsp$par))
env$N

dThing_dpars <- numderiv( log_final_abund, fit_hsp$par)
Vpar <- solve( fit_hsp$hessian) # the NEG() in optim() has already switched the sign
V_Thing <- dThing_dpars %*% (Vpar %*% dThing_dpars)
SE_log_final_abund <- sqrt( V_Thing)
CV_final_abund <- SE_log_final_abund
log_final_abund( fit_hsp$par) # reset the pop dyn in env to best-fitted


# ... and you can vectorize the (co)variance calcs, trivially, to give you covariance matrix

# More fancy diagnostics: compare E[POP] & Obs[POP] by some split of categories

## Here's another way to get the result, via a GLM
# Only works for exponential-in/decrease model
# Personally i prefer the "explicit lglk" approach, because you can put whatever you want
# into the pop dyn, rather than being constrained to the world of GLMs, and (ii) you get the pop dyn
# info immediately, But, sometimes a
# GLM can help you check (even if you have to approximate the model to make it GLMmable) and
# some audiences (NOT all) may feel more comfy with seeing GLM results--- at least til they are used to CKMR

# NB data is already reorgainzed into 'data.frame', during 'boring_data_prep_mammal'
# following was Eric A's suggestion. It is true that data.frame is easier to look at...
# ... but I was writing all this in C (or TMB...) like I would "for real", I'd stick with arrays

print( head( env$MHSP_df))

gfit <- glm(
    n_MHSP ~ B2 + dB + offset( log( n_comp_MHSP)),
    family=poisson( link=log),
    data=env$MHSP_df)

coef( gfit) # Should be MINUS the survival and RoI ests of the lglk version--- smaller N means *bigger*
fit_hsp$par
exp( fit_hsp$par[3]) # NB parametrization...

# If you were using a GLM "for real" yourself rather than the lglk, then YOU have to reconstruct what
# the actual N's were from the GLM coefficients...



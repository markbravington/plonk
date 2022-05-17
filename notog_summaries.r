## Stuff used for CAPAM 21 talk
## To be run only after fit_notog21.r

YEARS <- dimseq( env$N_ya)[[1]]
AFEC <- dimseq( env$N_ya)[[2]]
Vpar <- solve( fit$hessian)

# Which age makes most babies?
with( env, colSums( N_ya) * fecata[4:12])
# ... age 7

# TRO in equiv 7yo
TRO7_est_fun <- function( pars) { lglk_fish( pars); with( env,
    N_ya %**% fecata[ AFEC] / fecata[SLICE=7]
  )}
dTRO7 <- numderiv( TRO7_est_fun, fit$par)
TRO7_est <- TRO7_est_fun( fit$par) # ensures env has best vals by doing this after numderiv
V7 <- dTRO7 %**% Vpar %**% t( dTRO7)
SE7 <- sqrt( diag( V7))

# True value
TRO7_tru <- with( samp_notog1@secret,
    N_sya[ SLICE='F', YEARS, 4:12] %**% fec_sa[ SLICE='F',AFEC] / fec_sa[ SLICE='F',SLICE=7])

par( cex=1.5)
plot( YEARS, TRO7_tru, type='l', col='grey', ylim=c( 0, 1.1 * max( c( TRO7_est, TRO7_tru))),
  xlab='', ylab='')
title( main='TRO equiv 7yo', cex.main=1.5, line=-2)
points( YEARS, TRO7_est, col='darkred', pch=16)
error_bars( YEARS, TRO7_est, 1.96 * SE7, col='darkred')
n_POP_B <- sumover( env$n_POP, cq( Yad, Aad))
text( YEARS+0.3, TRO7_est, labels=unclass( n_POP_B[ YEARS]))

fec_tru <- with( samp_notog1@secret, fec_sa[ SLICE='F',AFEC] / fec_sa[ SLICE='F',SLICE=7])
fec_est_fun <- function( par) { lglk_fish( par); with( env,
    fecata[ AFEC] /fecata[SLICE=7]
  )}
dfec <- numderiv( fec_est_fun, fit$par)
fec_est <- fec_est_fun( fit$par)

Vfec <- dfec %**% Vpar %**% t( dfec)
SEfec <- sqrt( diag( Vfec))

plot( AFEC, fec_tru, type='l', col='grey', ylim=c( 0, 1.1 * max( c( fec_est + 2*SEfec, fec_tru))),
   xlab='', ylab='')
title(     main='Fec at age, rel to age 7', line= -2, cex.main=1.5)
points( AFEC, fec_est, col='blue', pch=16)
error_bars( AFEC, fec_est, 1.96 * SEfec, col='blue')

M_fun <- function( par) { lglk_fish( par); with( env,
    M
  )}

dM <- numderiv( M_fun, fit$par)
VM <- dM %**% Vpar %**% t( dM)
M_est <- M_fun( fit$par)

dSPRR <- numderiv( SPRR_fun, fit$par)
VSPRR <- dSPRR %**% Vpar %**% t( dSPRR)
SPRR_est <- SPRR_fun( fit$par)

# True SPRR: ensure same ranges as SPRR_est
SPRR_tru <- local({
  extract.named( SPRR_fun( fit$par, TRUE)) # create lastNyears and relevA
  SPRR <- with( samp_notog1@secret, {
    Zav <- colMeans( Z_sya[ SLICE='F', lastNyears, relevA])
    Mav <- M_sa[ SLICE='F', relevA]
    Na_now <- exp( -cumsum( c( 0, Zav)))
    Na_virgin <- exp( -cumsum( c( 0, Mav)))
    relevFec <- fec_sa[ SLICE='F', c( relevA, tail( relevA, 1)+1)]
    ROlife_virgin <- Na_virgin %*% relevFec
    ROlife_now <- Na_now %*% relevFec
    SPRR <- ROlife_now / ROlife_virgin
  return( c( SPRR))
  })
return( SPRR)
})


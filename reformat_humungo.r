# Loads the old Canada object and tweaks it...
# ... to match 2021 format

samp_humungo_A <- local({
  dataname <- load( 'd:/r2.0/canada/easywhale/nonUsamp_humungo1.rda') 
stopifnot( length( dataname)==1)
  get( dataname)
})

samp_humungo_A <- within( samp_humungo_A, {
  isamp <- 1:nrow( Samps)
  Samps$poss_par <- isamp %in% adid
  Samps$poss_off <- isamp %in% juvid
  MHSPs <- HSPs # ..I hope!
  rm( isamp, adid, juvid, POPs, HSPs)
})

save( samp_humungo_A, file='samp_humungo_A_raw.rda')


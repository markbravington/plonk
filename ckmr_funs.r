"numderiv" <-
function(f,x0,eps=0.0001, TWICE.=TRUE, param.name=NULL, ..., SIMPLIFY=TRUE) {
  if( is.null( param.name))
    ff <- function( x, ...) f( x, ...)
  else
    ff <- function( x, ...) {
      ll <- c( list( x), list( ...))
      names( ll)[1] <- param.name
      do.call( 'f', ll)
    }

  f0 <- ff(x0, ...)
  n <- length( x0)
  m <- matrix( 0, length(f0), n)
  for( i in 1:n) {
    this.eps <- eps * if( x0[ i]==0) 1 else x0[ i]
    m[,i] <- ( ff( x0+this.eps * (1:n==i), ...) - f0) / this.eps }
  if( !is.null( dim( f0)))
    dim( m) <- c( dim( f0), n)
  if( TWICE.) {
    mc <- match.call()
    mc$eps <- -eps
    mc$TWICE. <- FALSE
    m <- 0.5*(m + eval( mc, sys.frame( sys.parent())))
  }
  
  if( any( dim( m)==length( m)) && SIMPLIFY)
    m <- c( m) # demote 1D arrays to vectors
  
return( m)
}



".First.task" <-
function(...) {
  library( offarray)
  library( kinsimmer)
  if( !nzchar( Sys.which( 'g++'))) {
    makepath <- Sys.which( 'make')
    rtoolspath <- dirname( dirname( dirname( makepath)))
    bits <- sub( '.*[^0-9]([0-9]+)$', '\\1', R.version$arch)
    if( bits != '64') {
      bits <- 32 # fuck em all they don't fucken make it EASY do they ???
    }
    compiler_path <- file.path( rtoolspath, 'mingw' %&% bits, 'bin')
    if( !file.exists( file.path( compiler_path, 'g++.exe'))) {
      warning( "Can't find g++ in '" %&% compiler_path %&% "'")
    } else {
      Sys.setenv( PATH=sprintf( '%s;%s', compiler_path, Sys.getenv( 'PATH')))
    }
  }

  # All monkeys to typewriters! All monkeys to typewriters! Start typing...
  if( .Platform$r_arch=='x64' && nzchar( bp64 <- Sys.getenv( 'BINPREF64'))) {
    Sys.setenv( BINPREF=bp64)
  }

  library( TMB)
}



"add_handed" <-
function( df, pRight) {
  extract.named( df)

  my_Mum <- match( Mum, Me)
  my_Dad <- match( Dad, Me)

  n <- length( Mum)
  chirogene <- matrix( 0L, n, 2)
  is_founder <- Mum=='founder' # and Dad will be too, by defn
  chirogene[ is_founder,] <- rsample( 2*sum(is_founder), 1:2, prob=c( pRight, 1-pRight), replace=TRUE)
  # Everyone who's not a founder, gets their genes from an ancestor with knowable genes

  # Which Mat & which Pat gene does everyone get?
  whicho_Mum <- rsample( n, 1:2, replace=TRUE)
  whicho_Dad <- rsample( n, 1:2, replace=TRUE)

  repeat{
    nogenes <- which( rowSums( chirogene)==0)
    if( !length( nogenes)){
  break
    }

    # nkp is nogeners whose parent's genes *are* both now known
    nkp <- nogenes[ (chirogene[ my_Mum[ nogenes],1]>0) & (chirogene[ my_Dad[ nogenes],1]>0)]
    if( !length( nkp)) {
stop( "Impasse...") # no progress possible!
    }
    chirogene[ nkp,1] <- chirogene[ cbind( my_Mum[ nkp], whicho_Mum[ nkp])]
    chirogene[ nkp,2] <- chirogene[ cbind( my_Dad[ nkp], whicho_Dad[ nkp])]
  }

  Hand <- ifelse( rowSums( chirogene==1)==2, 'R', 'L')
  df$chirogene <- chirogene
  df$Hand <- Hand
return( df)
}



"ALT_generic_lglk_MOP_fuzzall_mammal" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  # ... that's all,  folks

  ## Population dynamics
  N <- autoloop( Y=years, {
    Nfad_y0 * exp( RoI * (Y-y0))
  })

  ## Ideal kinprob
  Pr_MOP_BYA <- autoloop( Bju=years, Ypar=Y_range, Apar=A_range, {
    Bpar <- Ypar - Apar
    ( Bju <= Ypar ) *          # was ad still alive at B[ju] ?
    ( Bju >= Amat + Bpar ) *   # was ad mature at B[ju] ?
    ( 1/N[ Bju])                  # ad's ERO at B[ju] if alive & mature, divided by TRO that year
  })

  ## Kinprob given OBSERVED covars, ie Fuzzage instead of true A
  
  ## Bayes' theorem to get Pr[ truA | Fuzzage, Y]
  Pr_AF_Y <- autoloop( A=A_range, F=F_range, Y=Y_range, {
    Pr_Fuzzage_Atru[ F, A] * 
    est_Pr_A_SampY[ A, Y]
  })
  
  Pr_F_Y <- sumover( Pr_AF_Y, 'A')
  Pr_A_FY <- autoloop( A=A_range, F=F_range, Y=Y_range, {
    Pr_AF_Y[ A, F, Y] / Pr_F_Y[ F, Y]
  })
    
  # "Mix" the ideal probs 

  if( 1 %in% calc_version){
    # Ideal autoloop version, but slow because the 6D array is HUGE!
    # Fancier sumover/autoloop combo function could fix this

		Pr_MOPAA_FYFY <- autoloop( A1=A_range, A2=A_range, F1=F_range, Y1=Y_range, F2=F_range, Y2=Y_range, {
			# True kinship (MO or OM) is unobservable. MO and OM are mutually exclusive, so we can add their probs...
			( Pr_MOP_BYA[ Y1-A1, Y2, A2] + Pr_MOP_BYA[ Y2-A2, Y1, A1]) *    
			 Pr_A_FY[ A1, F1, Y1] *
			 Pr_A_FY[ A2, F2, Y2]
		})
		
		Pr_MOP_FYFY <- sumover( Pr_MOPAA_FYFY, cq( A1, A2))
  }    
  
  # Faster alternative I hope. Keep the FOR-loops as small as possible...
  if( 2 %in% calc_version)evalq({
  	Pr_MOP_FYFY <- autoloop( F1=F_range, Y1=Y_range, F2=F_range, Y2=Y_range, 0) # store results
  	
		for( A1 in A_range) for( A2 in A_range) {
			this_Pr_MOP_FYFY <- autoloop( F1=F_range, Y1=Y_range, F2=F_range, Y2=Y_range, {
				( Pr_MOP_BYA[ Y1-A1, Y2, A2] + Pr_MOP_BYA[ Y2-A2, Y1, A1]) *    
				 Pr_A_FY[ A1, F1, Y1] *
				 Pr_A_FY[ A2, F2, Y2]
				})
			Pr_MOP_FYFY <- Pr_MOP_FYFY + this_Pr_MOP_FYFY
		}
  })

  if( 3 %in% calc_version)evalq({
  	Pr_MOP_FYFY <- autoloop( F1=F_range, Y1=Y_range, F2=F_range, Y2=Y_range, 0) # store results
  	
		for( F1 in F_range) for( F2 in F_range) {
			this_Pr_MOPAA_FYFY <- autoloop( A1=A_range, Y1=Y_range, A2=A_range, Y2=Y_range, {
				( Pr_MOP_BYA[ Y1-A1, Y2, A2] + Pr_MOP_BYA[ Y2-A2, Y1, A1]) *    
				 Pr_A_FY[ A1, F1, Y1] *
				 Pr_A_FY[ A2, F2, Y2]
				})
		  Pr_MOP_FYFY[ SLICE=F1,,SLICE=F2,] <- sumover( this_Pr_MOPAA_FYFY, cq( A1, A2))
		}
  })


  # Check: dimseq( Pr_MOP_FYFY) vs dimseq( n_MOP)
  Pr_MOP <- Pr_MOP_FYFY #  keep same notation as other examples

  # Useful to keep some stuff after function exits. Trust me, this works...
  list2env( mget( cq( Nfad_y0, RoI, N, Pr_MOP)), envir=environment( sys.function()))

  if( !all( is.finite( Pr_MOP) & (Pr_MOP >= 0) & (Pr_MOP <= 1))) {
return( -1e100) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  lglk <- sum( dpois( n_MOP, lambda=n_comp_MOP * Pr_MOP, log=TRUE))

return( lglk)
}



"boring_data_prep_delfi_A" <-
function( sampo, prev_env=parent.env()){
## Given "CKMR sample data", probably from 'prep_from_sim2'...
## Make arrays with n_comps and n_kin, per kinship and per covars
## Stick them into child env of prev_env--- which means you
## can subsequently associate them with a lglk function
## so it will "know" about its specific data

## usage:
# lglk_mydata <- generic_lglk_for_this_kinda_CKMR_data  # a function
# env <- boring_data_prep_mydata( my_samps, prev_env=environment( lglk_mydata))
# environment( lglk_mydata) <- env
# lglk_mydata( params) # lglk_mydata can now refer internally to 'n_MOP' etc

## Warning: this is boring. Did the name not tip you off? You almost
## certainly do NOT need to understand what's in it. And there's esoteRica
## which I am NOT going to explain. So you should probably stop reading RIGHT NOW...

## NB sampling has already ensured that only Females are taken--- so all POPs are actually MOPs
## and all Samps are Female


  extract.named( sampo) # Samps, POPs, HSPs, ...
  extract.named( Samps)
  extract.named( sampo@public) # Amat
  AMAX <- lastel( C_sya, 3)
  # Year of birth
  B <- Y - A
  
  adid <- which( poss_par)
  juvid <- which( poss_off)

  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number

  # Package up stuff, to be used as environment for lglk function
  y0 <- min( B[ juvid]) # SHOULDN'T really be data-driven
  years <- min( B[ juvid]) %upto% max( Y[ adid])
  envo <- list2env( mget( cq( rPOPs, B, Y, A, Amat, y0, years, juvid, adid)), parent=prev_env)

  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Bju_range <- min( B[ juvid]) %upto% max( B[ juvid])
  Yad_range <- min( Y[ adid]) %upto% max( Y[ adid])
  Aad_range <- min( A[ adid]) %upto%  max( A[ adid])

  # m_... is samp size
  m_ad_YA <- offarray( table( Y=Y[ adid], A=A[ adid]))
  m_ju_B <- offarray( table( B=B[ juvid]))

  # Number of comparisons: product of sample sizes by category
  # See ?offarray::noloop  or code of lglk_aggregate() below
  n_comp_MOP <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=Aad_range, {
      Bad <- Yad - Aad
      # Only do comps where ju is born after adult
      # ... which also avoids double-counting and self-comparisons
      m_ju_B[ Bju] * m_ad_YA[ Yad, Aad] * (Bju > Bad)
    })

  n_MOP <- offarray( table(
      Bju=B[ POPs[,2]],
      Yad=Y[ POPs[,1]],
      Aad=A[ POPs[,1]]),
      template=n_comp_MOP) # template ensures full ranges used, even if no entries for some values

  # That's what's needed for fitting...
  # ... next only for "pedagogical" purposes (data inspection)
  POP_df <- boring_dfize( n_MOP, n_comp_MOP)

  # copy useful vars into envo... R magic, just trust me on this one ;)
  list2env( mget( cq( n_comp_MOP, n_MOP,
      Bju_range, Yad_range, Aad_range, # in sample
      AMAX,
      POP_df)),
      envo)

return( envo)


}



"boring_data_prep_delfi_B" <-
function( sampo, prev_env=parent.env()){
## Given "CKMR sample data", probably from 'prep_from_sim2'...
## Make arrays with n_comps and n_kin, per kinship and per covars
## Stick them into child env of prev_env--- which means you
## can subsequently associate them with a lglk function
## so it will "know" about its specific data

## usage:
# lglk_mydata <- generic_lglk_for_this_kinda_CKMR_data  # a function
# env <- boring_data_prep_mydata( my_samps, prev_env=environment( lglk_mydata))
# environment( lglk_mydata) <- env
# lglk_mydata( params) # lglk_mydata can now refer internally to 'n_MOP' etc

## Warning: this is boring. Did the name not tip you off? You almost
## certainly do NOT need to understand what's in it. And there's esoteRica
## which I am NOT going to explain. So you should probably stop reading RIGHT NOW...

# This one has a "full" version of the data, plus a "_noHands" version which *should* give
# biased estimates under a naive model (or is inestimable under the right model!)

## NB sampling has already ensured that only Females are taken--- so all POPs are actually MOPs
## and all Samps are Female


  extract.named( sampo) # Samps, POPs, HSPs, ...
  extract.named( Samps)
  extract.named( sampo@public) # Amat
  AMAX <- lastel( C_sya, 3)
  # Year of birth
  B <- Y - A

  adid <- which( poss_par)
  juvid <- which( poss_off)

  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number

  # Package up stuff, to be used as environment for lglk function
  y0 <- min( B[ juvid]) # SHOULDN'T really be data-driven
  years <- min( B[ juvid]) %upto% max( Y[ adid])
  envo <- list2env( mget( cq( rPOPs, B, Y, A, Amat, y0, years, juvid, adid)), parent=prev_env)

  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Bju_range <- min( B[ juvid]) %upto% max( B[ juvid])
  Yad_range <- min( Y[ adid]) %upto% max( Y[ adid])
  Aad_range <- min( A[ adid]) %upto%  max( A[ adid])

  # m_... is samp size
  m_ad_YAH <- offarray( table( Y=Y[ adid], A=A[ adid], H=Hand[ adid]))
  m_ju_BH <- offarray( table( B=B[ juvid], H=Hand[ juvid]))
  HANDS <- dimseq( m_ju_BH)$H

  # Number of comparisons: product of sample sizes by category
  # See ?offarray::noloop  or code of lglk_aggregate() below
  n_comp_MOP <- autoloop( Bju=Bju_range, Hju=HANDS, Yad=Yad_range, Aad=Aad_range, Had=HANDS, {
      Bad <- Yad - Aad
      # Only do comps where ju is born after adult
      # ... which also avoids double-counting and self-comparisons
      m_ju_BH[ Bju, Hju] * m_ad_YAH[ Yad, Aad, Had] * (Bju > Bad)
    })

  n_MOP <- offarray( table(
      Bju=B[ POPs[,2]], Hju=Hand[ POPs[,2]],
      Yad=Y[ POPs[,1]], Aad=A[ POPs[,1]], Had=Hand[ POPs[,1]]),
      template=n_comp_MOP) # template ensures full ranges used, even if no entries for some values


  # Hands-free versions
  m_ad_YA <- sumover( m_ad_YAH, 'H')
  m_ju_B <- sumover( m_ju_BH, 'H')
  n_comp_MOP_noHands <- sumover( n_comp_MOP, cq( Had, Hju))
  n_MOP_noHands <- sumover( n_MOP, cq( Had, Hju))

  # That's what's needed for fitting...
  # ... next only for "pedagogical" purposes (data inspection)
  MOP_df <- boring_dfize( n_MOP, n_comp_MOP) 
  MOP_noHands_df <- boring_dfize( n_MOP_noHands, n_comp_MOP_noHands) 

  # copy useful vars into envo... R magic, just trust me on this one ;)
  list2env( mget( cq( n_comp_MOP, n_MOP, n_MOP_noHands, n_comp_MOP_noHands,
      Bju_range, Yad_range, Aad_range, # in sample
      AMAX, HANDS,
      MOP_df, MOP_noHands_df)), envo)

return( envo)


}



"boring_data_prep_delfi_C" <-
function( sampo, prev_env=parent.env()){
## Given "CKMR sample data", probably from 'prep_from_sim2'...
## Make arrays with n_comps and n_kin, per kinship and per covars
## Stick them into child env of prev_env--- which means you
## can subsequently associate them with a lglk function
## so it will "know" about its specific data

## usage:
# lglk_mydata <- generic_lglk_for_this_kinda_CKMR_data  # a function
# env <- boring_data_prep_mydata( my_samps, prev_env=environment( lglk_mydata))
# environment( lglk_mydata) <- env
# lglk_mydata( params) # lglk_mydata can now refer internally to 'n_MOP' etc

## Warning: this is boring. Did the name not tip you off? You almost
## certainly do NOT need to understand what's in it. And there's esoteRica
## which I am NOT going to explain. So you should probably stop reading RIGHT NOW...

## NB sampling has already ensured that only Females are taken--- so all POPs are actually MOPs
## and all Samps are Female

  require( deconvodisc) # might as well stop early if it's not present
  extract.named( sampo) # Samps, POPs, HSPs, ...
  extract.named( Samps)
  extract.named( sampo@public) # Amat, fuzzer
  AMAX <- lastel( C_sya, 3)

  adid <- which( poss_par)
  juvid <- which( poss_off)

  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number

  ## Package up stuff, to go into envir of lglk function
  B <- Y - Jage # will be NA for adults, cos only Fuzzage avail
  y0 <- min( B[ juvid]) # SHOULDN'T really be data-driven
  years <- min( B[ juvid]) %upto% max( Y[ adid])
  envo <- list2env( mget( cq( rPOPs, B, Y, Fuzzage, Amat, y0, years, juvid, adid)), parent=prev_env)

  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Bju_range <- min( B[ juvid]) %upto% max( B[ juvid])
  Yad_range <- min( Y[ adid]) %upto% max( Y[ adid])
  Aad_range <- dimseq( fuzzer)$Atru
  Fuzzad_range <- min( Fuzzage[ adid]) %upto%  max( Fuzzage[ adid])

  # m_... is samp size
  m_ad_YFuzz <- offarray( table( Y=Y[ adid], Fuzz=Fuzzage[ adid]))
  m_ju_B <- offarray( table( B=B[ juvid]))

  # Number of comparisons: product of sample sizes by category
  # See ?offarray::noloop  or code of lglk_aggregate() below
  n_comp_MOP <- autoloop( Bju=Bju_range, Yad=Yad_range, Fuzzad=Fuzzad_range, {
      # Do _all_ comps --- don't worry about 2nd-guessing adult maturity
      # which I did at first, but got it wrong!
      # Pr_MOP should sort out impossibles
      m_ju_B[ Bju] * m_ad_YFuzz[ Yad, Fuzzad]
    })

  n_MOP <- offarray( table(
      Bju=B[ POPs[,2]],
      Yad=Y[ POPs[,1]],
      Fuzzad=Fuzzage[ POPs[,1]]),
      template=n_comp_MOP) # template ensures full ranges used, even if no entries for some values

  # Now estimate *true* age compo in *samples*, by deconvolution of fuzzage compo
  samp_Fuzzad <- offarray( table(
    Fuzzage=Fuzzage[ adid],
    Y=Y[ adid]
  ))

  # fuzzer has more age-classes than seen in actual data, so trim it...
  # ... 1 SD in from max Fuzzage. This is AD HOC but shouldn't matter; no good exact solution.
  # LONG BORING SNIPPET !!!
  # NB we are deliberately using the _fuzzed_ age to slice in the _true_ age dimension
  # which is Weird. Point is that the SD should not vary that much in the older years
  #  # Hope that's clear :)
  MFA <- tail( Fuzzad_range, 1) # Max Fuzz Age seen

  # Want SD of fuzzed age at that (true) age. Use complete span of allowed Fuzzages for this
  Fuzzage_range <- dimseq( fuzzer)[[1]] # actually called 'Aconv' for some reason
  fuzz_MFA <- fuzzer[ Fuzzage_range, SLICE=MFA]
  EEE0 <- sum( fuzz_MFA) # should be 1 AFAICS
  EEE1 <- (fuzz_MFA %*% Fuzzage_range) / EEE0
  EEE2 <- (fuzz_MFA %*% sqr( Fuzzage_range)) / EEE0
  SD1 <- sqrt( EEE2 - sqr( EEE1))
  likely_Amax <- c( floor( MFA - SD1)) # c() to de-matrix
  trim_Aad_range <- Aad_range %such.that% (. < likely_Amax)
  
  # WARNING--  I found this unstable. Changing the "<" to "<=" above, gives much wobblier estimates...
  Pr_Fuzzage_Atru <- fuzzer[ Fuzzad_range, trim_Aad_range]
  names( dimnames( Pr_Fuzzage_Atru)) <- cq( Fuzzage, A)
  # Normalize columns to sum to 1
  sumbo <- sumover( Pr_Fuzzage_Atru, 'Fuzzage')
  Pr_Fuzzage_Atru <- autoloop( Fuzzage=Fuzzad_range, A=trim_Aad_range, {
    Pr_Fuzzage_Atru[ Fuzzage, A] / sumbo[ A]
  })

  # deconvodisc() expects normal R matrix, hence unclass()
  suppressWarnings({ # prolly shouldn't in general; here nlminb gives harmless warning...
    # ... don't wanna frighten people
    temp <- deconvodisc::deconvodisc( samp_Fuzzad, unclass( Pr_Fuzzage_Atru))
  })
  # ... ie est distros, plus covar mat of each one (by year). Need the former, with correct dims:
  est_Pr_A_SampY <- offarray( temp$pr_truhat,
      first=c( A= head( trim_Aad_range, 1), Y=firstel( samp_Fuzzad, 'Y')),
      last=c( A= tail( trim_Aad_range, 1), Y=lastel( samp_Fuzzad, 'Y')))

  # That's what's needed for fitting...
  # ... next only for "pedagogical" purposes (data inspection)
  MOP_df <- boring_dfize( n_MOP, n_comp_MOP) 

  # copy useful vars into envo... R magic, just trust me on this one ;)
  likely_Aad_range <- trim_Aad_range # name used in lglk
  list2env( mget( cq( 
      n_comp_MOP, n_MOP,
      Bju_range, Yad_range, Fuzzad_range, likely_Aad_range,
      est_Pr_A_SampY,
      Pr_Fuzzage_Atru,
      samp_Fuzzad,
      AMAX,
      MOP_df)),
    envo)

return( envo)


}



"boring_data_prep_delfi_D" <-
function( sampo, prev_env=parent.env()){
## Given "CKMR sample data", probably from 'prep_from_sim2'...
## Make arrays with n_comps and n_kin, per kinship and per covars
## Stick them into child env of prev_env--- which means you
## can subsequently associate them with a lglk function
## so it will "know" about its specific data

## usage:
# lglk_mydata <- generic_lglk_for_this_kinda_CKMR_data  # a function
# env <- boring_data_prep_mydata( my_samps, prev_env=environment( lglk_mydata))
# environment( lglk_mydata) <- env
# lglk_mydata( params) # lglk_mydata can now refer internally to 'n_POP' etc

## Warning: this is boring. Did the name not tip you off? You almost
## certainly do NOT need to understand what's in it. And there's esoteRica
## which I am NOT going to explain. So you should probably stop reading RIGHT NOW...

## NB sampling has already ensured that only Females are taken--- so all POPs are actually MOPs
## and all Samps are Female (but PHSPs are still noted, tho not used in this one)

  require( deconvodisc) # might as well stop early if it's not present
  extract.named( sampo) # Samps, POPs, HSPs, ...
  extract.named( Samps) # Y, Fuzzage
  extract.named( sampo@public) # Amat, fuzzer
  AMAX <- lastel( C_sya, 3)
  PLAUS_AGE_TRIMMER <- 0.02 # discard (older) ages whose rev cumul prob falls below this

  # This time *everyone* is a potential parent *and* a potential offspring...
  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number

  ## Package up stuff, to go into envir of lglk function
  y0 <- min(Y)-Amat       # Pretty arbitrary! Interested only from here on
  ylast <- max(Y)
  years <- y0 %upto% max( Y) # for pop dyn model. This could go waaay back in time than, but we will truncate
  samp_Y_range <- min( Y) %upto% max( Y) 
  
  envo <- list2env( mget( cq( rPOPs, Y, Fuzzage, Amat, y0, ylast, years, samp_Y_range)), parent=prev_env)

  # Now estimate *true* age compo in *samples*, by deconvolution of fuzzage compo
  samp_FuzzY <- offarray( table(
    Fuzzage=Fuzzage,
    Y=Y
  ))
  
  # fuzzer has more age-classes than seen in actual data, so trim it...
  # ... 1 SD in from max Fuzzage. This is AD HOC but shouldn't matter; no good exact solution.
  # LONG BORING SNIPPET !!!
  # NB we are deliberately using the _fuzzed_ age to slice in the _true_ age dimension
  # which is Weird. Point is that the SD should not vary that much in the older years
  #  # Hope that's clear :)
  MFA <- max( Fuzzage) # Max Fuzz Age seen
  samp_Fuzz_range <- 1 %upto% MFA

  # Want SD of fuzzed age at that (true) age. Use complete span of allowed Fuzzages for this
  Fuzzage_range <- dimseq( fuzzer)[[1]] # actually called 'Aconv' for some reason
  fuzz_MFA <- fuzzer[ Fuzzage_range, SLICE=MFA]
  EEE0 <- sum( fuzz_MFA) # should be 1 AFAICS
  EEE1 <- (fuzz_MFA %*% Fuzzage_range) / EEE0
  EEE2 <- (fuzz_MFA %*% sqr( Fuzzage_range)) / EEE0
  SD1 <- sqrt( EEE2 - sqr( EEE1))
  likely_Amax <- c( floor( MFA - SD1)) # c() to de-matrix
  trim_A_range <- 1 %upto% likely_Amax
  
  Pr_Fuzzage_Atru <- fuzzer[ samp_Fuzz_range, trim_A_range]
  names( dimnames( Pr_Fuzzage_Atru)) <- cq( Fuzzage, A)
  # Normalize columns to sum to 1
  sumbo <- sumover( Pr_Fuzzage_Atru, 'Fuzzage')
  Pr_Fuzzage_Atru <- autoloop( Fuzzage=samp_Fuzz_range, A=trim_A_range, {
    Pr_Fuzzage_Atru[ Fuzzage, A] / sumbo[ A]
  })

  # deconvodisc() expects normal R matrix, hence unclass()
  suppressWarnings({ # prolly shouldn't in general; here nlminb gives harmless warning...
    # ... don't wanna frighten people
    temp <- deconvodisc::deconvodisc( samp_FuzzY, unclass( Pr_Fuzzage_Atru))
  })
  # ... ie est distros, plus covar mat of each one (by year). Need the former, with correct dims:
  est_Pr_A_SampY <- offarray( temp$pr_truhat,
      first=c( A= head( trim_A_range, 1), Y=head( samp_Y_range,1)),
      last=c( A= tail( trim_A_range, 1), Y=tail(samp_Y_range,1)))


  ## Stuff for aggregated version:
  
  ## Only include comps where AT LEAST ONE of the pair is very likely born afterat y0
  # MASSIVE PITA
  # Here's the code from the lglk...
	Pr_AF_Y <- autoloop( A=trim_A_range, F=samp_Fuzz_range, Y=samp_Y_range, {
		Pr_Fuzzage_Atru[ F, A] * 
		est_Pr_A_SampY[ A, Y]
	})

	Pr_F_Y <- sumover( Pr_AF_Y, 'A')
	Pr_A_FY <- autoloop( A=trim_A_range, F=samp_Fuzz_range, Y=samp_Y_range, {
		Pr_AF_Y[ A, F, Y] / Pr_F_Y[ F, Y]
	})

  ## Now we need min *plaus* B for each F & Y
  # Find cumul probs of A given F, Y; then find the one that's within PLAUS_AGE_TRIMMER of 1
  # Probably autoloopable, but tricky and unclear...
  min_plaus_B_FY <- autoloop( F=samp_Fuzz_range, Y=samp_Y_range, 0) # easier than calling offarray() !
  evalq({ # for debugging speed
		for( iF in samp_Fuzz_range){
			for( iY in samp_Y_range) {
				cumPr_A_FY <- cumsum( unclass( Pr_A_FY[,SLICE=iF,SLICE=iY]))
				max_plaus_A_FY <- 1 + findInterval( 1-PLAUS_AGE_TRIMMER, cumPr_A_FY)
				min_plaus_B_FY[ iF, iY] <- iY - max_plaus_A_FY
			}
		}
  })
  
  # m_... is samp size
  m_FY <- samp_FuzzY

  # This lglk needs n_MOP_FYFY
  # Number of comparisons: product of sample sizes by category
  
  # POPs will be organized so that Y1 < Y2, or F1 <= F2 if Y1=Y2
  # When Y1==Y2 and F1==F2, do 1/2 the number of comps (really m*(m-1)/2 but...)
  # Thus, only do comparisons that correspond
  # Also will discard POPs where O would plausibly be prior to y0
  n_comp_MOP <- autoloop( F1=samp_Fuzz_range, Y1=samp_Y_range, F2=samp_Fuzz_range, Y2=samp_Y_range, {
      ((min_plaus_B_FY[ F1, Y1] >= y0) | (min_plaus_B_FY[ F2, Y2] >= y0)) * 
      m_FY[ F1, Y1] * m_FY[ F2, Y2] *  
      ( ((Y1<Y2) | ((Y1==Y2) & (F1<F2))) + 0.5 * ((Y1==Y2) & (F1==F2)) )
    })

  # Reorder the POPs so 
	F1 <- Fuzzage[ POPs[,1]]
	Y1 <- Y[ POPs[,1]]
	F2 <- Fuzzage[ POPs[,2]]
	Y2 <- Y[ POPs[,2]]
  swappo <- (Y1>Y2) | ((Y1==Y2) & (F1>F2))
  POPs[ swappo, 1:2] <- POPs[ swappo, 2:1]
 
	# check that POPs have been reorganized OK!
	F1 <- Fuzzage[ POPs[,1]]
	Y1 <- Y[ POPs[,1]]
	F2 <- Fuzzage[ POPs[,2]]
	Y2 <- Y[ POPs[,2]]
stopifnot( all( ((Y1<Y2) | ((Y1==Y2) & (F1<F2))) | ((Y1==Y2) & (F1==F2)) ))


  POPs <- POPs[( (min_plaus_B_FY[ MATSUB=cbind( F1, Y1)] >= y0) | (min_plaus_B_FY[ MATSUB=cbind( F2, Y2)] >= y0) ),]

  n_MOP <- offarray( table(
      F1=Fuzzage[ POPs[,1]],
      Y1=Y[ POPs[,1]],
      F2=Fuzzage[ POPs[,2]],
      Y2=Y[ POPs[,2]]),
      template=n_comp_MOP) # template ensures full ranges used, even if no entries for some values

  if( TRUE || paranoid) {
stopifnot( all( n_MOP[ VECSUB=which( c( n_comp_MOP)==0)] == 0) )
  }

  # That's what's needed for fitting...
  # ... next only for "pedagogical" purposes (data inspection)
  MOP_df <- boring_dfize( n_MOP, n_comp_MOP) 

  # copy useful vars into envo... R magic, just trust me on this one ;)
  A_range <- trim_A_range # for lglk
  F_range <- samp_Fuzz_range
  Y_range <- samp_Y_range
  list2env( mget( cq( 
      n_comp_MOP, n_MOP,
      Y_range,F_range, A_range,
      est_Pr_A_SampY,
      Pr_Fuzzage_Atru,
      AMAX,
      MOP_df)),
    envo)

return( envo)


}



"boring_data_prep_fish" <-
function( sampo, prev_env=parent.env()){
  extract.named( sampo)
  extract.named( Samps)
  extract.named( sampo@public) # Amat, catches
  # Year of birth
  B <- Y - A

  # POPs (and HSPs) are supposed to be in birth-order
  # ... with same-cohort HSPs already excluded
  # Obvs couldn't do things exactly this way for real data with age uncertainty
  # Could adjust order etc here, but earlier data prep is supposed to ensure it, so it's a check on misprep...
stopifnot(
    all( B[ POPs[,1]] < B[ POPs[,2]]),
    all( B[ HSPs[,1]] < B[ HSPs[,2]])
  )

  # For now, keep same-cht HSPs--- zap later

  # ONLY use POPS where juve is born *before* year-of-adult-capture
  # (to avoid bias when adults caught within-spawning-season don't get a fair go that year)
  # Strictly > (not >=)
  POPs <- POPs[ Y[ POPs[,1]] > B[ POPs[,2]],]

  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number

  # Package up stuff, to be used as environment for lglk function
  y0 <- min( B[ juvid]) # SHOULDN'T really be data-driven
  years <- min( B[ juvid]) %upto% max( Y[ adid])
  envo <- list2env( mget( cq( rPOPs, B, Y, A, Amat, y0, years, juvid, adid)), parent=prev_env)

  ## Stuff for aggregated version:
  # NB things will be evaluate over entire range, even if gaps
  # Dodge-able at C level (ie gaps could be handled), but maybe not in R
  Bju_range <- min( B[ juvid]) %upto% max( B[ juvid])
  Yad_range <- min( Y[ adid]) %upto% max( Y[ adid])
  Aad_range <- min( A[ adid]) %upto% max( A[ adid])

  # m_... is samp size
  m_ad_YA <- offarray( table( Y[ adid], A[ adid]))
  m_ju_B <- offarray( table( B[ juvid]))
  m_ad_Y <- offarray( table( Y[ adid])) # for noage

  # Number of comparisons: product of sample sizes by category
  # See ?offarray::autoloop  or code of lglk_aggregate() below
  n_comp_POP <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=Aad_range, {
      Bad <- Yad - Aad
      # Only do comps where ju is born after adult
      # ... which also avoids double-counting and self-comparisons
      # Also disallow year-of-adult-cap
      m_ju_B[ Bju] * m_ad_YA[ Yad, Aad] * (Bju > Bad) * (Bju < Yad)
    })


  # Now tot up number of POPs seen, in the same way
  n_POP <- offarray( table( Bju=B[ POPs[,2]], Yad=Y[ POPs[,1]], Aad=A[ POPs[,1]]),
      template=n_comp_POP)


  # And for HSPs... guaranteed in birth-order by row, thx2 prepare_from_sim()
  n_comp_HSP <- autoloop( B1=Bju_range, B2=Bju_range, {
      # NB *exclude* double-count and same-cohort
      m_ju_B[ B1] * m_ju_B[ B2] * (B2>B1)
    })

  n_HSP <- offarray( table( B1=B[ HSPs[,1]], B2=B[ HSPs[,2]]),
      template=n_comp_HSP)

  n_comp_POP_noage <- autoloop( Bju=Bju_range, Yad=Yad_range, {
      # No point if ju born after adult was (lethally) sampled
      m_ju_B[ Bju] * m_ad_Y[ Yad] * (Bju <= Yad)
    })


  # Now tot up number of POPs seen, in the same way
  n_POP_noage <- offarray( table( Bju=B[ POPs[,2]], Yad=Y[ POPs[,1]]),
      template=n_comp_POP_noage)

  # Age-structured stuff
  # poss discrep between "adults" in the test/sample and Amat
  AMAX <- max( Aad_range)

  min_Bju <- min( Bju_range)
  max_Bju <- max( Bju_range)

  # copy useful vars into envo... R magic
  weight <- attr( sampo, 'public')$wt_a
  Cfem_ya <- attr( sampo, 'public')$C_sya[ SLICE='F', min_Bju %upto% max_Bju, 1:AMAX]
  AminC <- min( which( colSums( Cfem_ya) > 0)) # start pop dyn from this age up
  ALL_AGES <- AminC %upto% AMAX
  Anonplus <- head( ALL_AGES, -1)
  ADAGES <- Amat %upto% AMAX
  Pr_FNeg_HSP <- 0 # usually ~0.15 in our real cases; here would need to "thin" true HSPs accordingly;
  # ... too lazy

  SIGNAL_BADFIT <- (-.Machine$double.xmax * 2^-32) # either optim or nlminb is sniffy about -Inf. The 2^blah is
    # ... to avoid getting toooo close to overflow in case of internal optim/nlminb calcs. Paranoid...

  list2env( mget( cq( n_comp_POP, n_POP, n_comp_HSP, n_HSP, n_POP_noage, n_comp_POP_noage,
      Bju_range, Yad_range, Aad_range,
      ALL_AGES, ADAGES, AMAX, Anonplus, AminC, min_Bju, max_Bju,
      weight, Cfem_ya,
      Pr_FNeg_HSP,
      SIGNAL_BADFIT
      )), envo)

  envo$last_params <- NA # for debugging
return( envo)
}



"boring_data_prep_humungo_A" <-
function( sampo, prev_env=parent.env()){
  extract.named( sampo) # Samps, POPs, MHSPs, ...
  extract.named( Samps) # Y, Fuzzage
  extract.named( sampo@public) # Amat, fuzzer
  AMAX <- lastel( C_sya, 3)

  juvid <- which( poss_off)
  
  B <- Y-A
	y0 <- min( B[ juvid]) # SHOULDN'T really be data-driven
	ylast <- max( B[ juvid]) 
	years <- y0 %upto% ylast
	Bju_range <- years
	
	envo <- list2env( mget( cq( B, Y, A, Amat, y0, ylast, years, juvid)), parent=prev_env)

  # And for MHSPs... guaranteed in birth-order by row, thx2 prepare_from_sim()
  m_ju_B <- offarray( table( B[ juvid]))
  n_comp_MHSP <- autoloop( B1=Bju_range, B2=Bju_range, {
      # NB *exclude* double-count and same-cohort
      m_ju_B[ B1] * m_ju_B[ B2] * (B2>B1)
    })

  # Make sure first-born is first
  swappo <- B[ MHSPs[,1]] > B[ MHSPs[,2]]
  MHSPs[ swappo,] <- MHSPs[ swappo,2:1]
  
  # Elim same-cohort MHSPs--- though mostly same-chters would be FSPs
  droppo <- B[ MHSPs[,1]] == B[ MHSPs[,2]]
  MHSPs <- MHSPs[ !droppo,]

  n_MHSP <- offarray( table( B1=B[ MHSPs[,1]], B2=B[ MHSPs[,2]]),
      template=n_comp_MHSP)

  # That's what's needed for fitting
  # For "pedagogical" purposes, convert those arrays into data.frames
  # this uses code from the base-R and the "markiverse", not the "tidyverse" !
  # ... and, I might add, it works just fine ;)

  MHSP_df <- as.data.frame( n_comp_MHSP, name_of_response='n_comp_MHSP')
  temp <- as.data.frame( n_MHSP, name_of_response='n_MHSP')
  MHSP_df <- cbind( MHSP_df, n_MHSP=temp$n_MHSP) # really ought to check that the index rows match...
  MHSP_df$dB <- with( MHSP_df, B2-B1)
  # Just the useful bits:
  MHSP_df <- MHSP_df %where% (n_comp_MHSP > 0)

  list2env( mget( cq( n_comp_MHSP, n_MHSP,
      Bju_range, # in sample
      AMAX, # why? who? wot?
      MHSP_df)), envo)

return( envo)
}



"boring_data_prep_rhombi_D" <-
function( sampo, prev_env=parent.env()){
  SEXES <- cq( F, M)
  envo <- list2env( mget( 'SEXES'), parent=prev_env)

  extract.named( sampo) # Samps, POPs, HSPs, ...
  extract.named( Samps)
  extract.named( sampo@public) # Amat
  
  adid <- which( poss_par)
  juvid <- which( poss_off)

  # Can't remember what this is for!
  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number
  
  # Not using offarray( table(...)) here becos it seems broken on 1D case!
  m_ad_S <- c( F=sum( Sex[ adid]=='F'), M=sum( Sex[ adid]=='M'))
  m_ju <- length( juvid) # or sum( is.na( Sex)), or...

  # Number of comparisons: product of sample sizes by category
  # See ?offarray::noloop  or code of lglk_aggregate() below
  n_comp_POP <- m_ju * m_ad_S

  
  # Following doesn't work becos offarray( table()) bug in 1D case...
#  n_POP <- offarray( table(
#      Sad=Sex[ POPs[,1]]),
#      template=n_comp_POP) # template ensures full ranges used, even if no entries for some values

  n_POP <- c( F=sum( Sex[ POPs[,1]]=='F'), M=sum( Sex[ POPs[,1]]=='M'))

  list2env( mget( cq( 
      n_comp_POP, n_POP)),
    envo)

return( envo)
}



"boring_data_prep_rhombi_W" <-
function( sampo, prev_env=parent.env()){
  SEXES <- cq( F, M)
  envo <- list2env( mget( 'SEXES'), parent=prev_env)

  extract.named( sampo) # Samps, POPs, HSPs, ...
  extract.named( Samps)
  extract.named( sampo@public) # Co_samp_m_UB
  Pr_Co_m_unbiased <- Co_samp_m_UB / sum( Co_samp_m_UB)
  Pr_Co_unbiased <- rbind( Co_levels==1, Pr_Co_m_unbiased)
  dimnames( Pr_Co_unbiased) <- list( S=SEXES, Co=NULL) # doesn't need to be offarray
  Pr_Co_unbiased@offset <- c( NA, NA)
  oldClass( Pr_Co_unbiased) <- 'offarray'
  
  adid <- which( poss_par)
  juvid <- which( poss_off)

  rPOPs <- pairid( POPs[,1], POPs[,2]) # combine into single real number
  
  m_ju <- length( juvid) # or sum( is.na( Sex)), or...
  m_ad_SCo <- offarray( table( S=Sex[ adid], Co=Co[ adid]))
  
  # Number of comparisons: product of sample sizes by category
  n_comp_POP <- m_ju * m_ad_SCo

  n_POP <- offarray( table(
      S=Sex[ POPs[,1]],
      Co=Co[ POPs[,1]]),
      template=n_comp_POP) # template ensures full ranges used, even if no entries for some values

  list2env( mget( cq( 
      Co_levels, n_Co_levels,
      Pr_Co_unbiased,
      n_comp_POP, n_POP)),
    envo)

return( envo)
}



"boring_dfize" <-
function( n_MOP, n_comp_MOP) {
  # Offarrays into data.frames
  # Uses code from base-R and the "markiverse", not the "tidyverse" !
  # ... and, I might add, it works just fine ;)
  MOP_df <- as.data.frame( n_comp_MOP, name_of_response='n_comp_MOP')
  temp <- as.data.frame( n_MOP, name_of_response='n_MOP')
  MOP_df <- cbind( MOP_df, n_MOP=temp$n_MOP) # really ought to check that the index rows match...
  MOP_df <- MOP_df %where% (n_comp_MOP > 0)
return( MOP_df)
}



"cheat_ruff_tru_pars" <-
function( simpop, bdp_env) {
  # bdp_env from boring_data_prep
  # defined in this crazy way to avoid having to preset envir of this fun

  funge <- function( sim) {
    # extract.named( make_sim_names()) # not needed AFAICS
    public <- attr( sim, 'public')
    secret <- attr( sim, 'secret')
    Amat <- public$Amat
    AMAX <- lastel( public$C_sya, 3)
    AGE <- 1:AMAX
    AminC <- min( which( colSums( public$C_sya[ SLICE='F',,]) > 0)) # start pop dyn from this age up

    # NB NB: *females* *only*
    # c() to strip offarray attrib
    fitto <- lm( log( c( secret$N_sya[ SLICE='F', SLICE=y0,])) ~ AGE, subset= AGE >= AminC)
    aslope_y0 <- abs( coef( fitto)[2])
    # 2021: logical subscript seems to kybosh the next
    # Ncad_y0 <- sum( secret$N_sya[ SLICE='F', SLICE=y0, AGE >= AminC])
    AGE_at_least_AminC <- AGE %such.that% (. >= AminC)
    Ncad_y0 <- sum( secret$N_sya[ SLICE='F', SLICE=y0, AGE_at_least_AminC])
    AGE_at_least_Amat <- AGE %such.that% (. >= Amat)
    M <- mean( secret$M[ SLICE='F', AGE_at_least_Amat]) # ish...
    fecata <- unclass( secret$fec_sa[ SLICE='F',]) # dropping
    fitto <- lm( log( fecata) ~ log( public$wt_a), subset=AGE >= Amat)
    wtfecpar <- coef( fitto)[2]

    pars <- mget( cq( Ncad_y0, aslope_y0, M, wtfecpar))
    pars <- sapply( pars, unname) # otherwise you get weird stuff tacked on
  return( log( pars))
  }
  environment( funge) <- bdp_env

  # make life easy when I'm trying to debug this...
  if( all( cq( debug, evaluator) %in% all.names( body( sys.function())[[2]]))) {
    mtrace( funge)
  }

return( funge( simpop))
}



"error_bars" <-
function( x, ybar, yadd, ysub=yadd, linmul=0.2, ...) {
  inches <- par( 'pin')[1] / diff( par('usr')[1:2])
  lin <- linmul * inches / length( x)
  suppressWarnings(  # about zero-length arrows, yawn
    arrows( x0=x, y0=ybar + yadd, x1=x, y1=ybar-ysub, code=3, angle=90, length=lin, ...)
  )
}



"fec_est_fun" <-
function( par) { lglk_fish( par); with( env,
    fecata[ AFEC] /fecata[SLICE=7]
  )}



"generic_lglk_cartoon" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfemale <- exp( params[ 1])
  Nmale <- exp( params[ 2])

  ## "Population dynamics" 
  Nad <- c( F=Nfemale, M=Nmale)

  ## Ideal kinprob
  Pr_POP_Sad <- 1 / Nad
  Pr_POP <- Pr_POP_Sad # for consistency with other generic lglks

  list2env( mget( cq( Nad, Pr_POP)), envir=environment( sys.function()))

  if( !all( is.finite( Pr_POP) & (Pr_POP >= 0) & (Pr_POP <= 1))) {
return( -1e100) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  lglk <- sum( dpois( n_POP, lambda=n_comp_POP * Pr_POP, log=TRUE))

return( lglk)
}



"generic_lglk_Covar" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfemale <- exp( params[ 1])
  Nmale <- exp( params[ 2])
  Co_expo <- c( F=0, M=params[3]) # hardwire no effect for Females 

  ## "Population dynamics" and "biology"
  Nad <- c( F=Nfemale, M=Nmale)
  fec_Co <- autoloop( S=SEXES, Co=Co_levels, {
    Co ^ Co_expo[S]
  })
  fec_Co_m <- unclass( fec_Co[ SLICE='M',])
  
  ## Compute TRO. We "know" true ppn by Co-level in adult Males thx2 unbiased sample...
  RO_SCo <- autoloop( S=SEXES, Co=Co_levels, {
    Nad[ S] * 
    fec_Co[ S, Co] * 
    Pr_Co_unbiased[ S, Co]
  })
  TRO <- sumover( RO_SCo, 'Co') # by Sex

  ## Ideal kinprob
  Pr_POP_SCo <- autoloop( S=SEXES, Co=Co_levels, {
    fec_Co[ S, Co] / TRO[ S]
  })
  
  ## Housekeeping...
  Pr_POP <- Pr_POP_SCo # for consistency with other generic lglks
  list2env( mget( cq( Nad, Pr_POP, Co_expo, fec_Co_m)), envir=environment( sys.function()))

  if( !all( is.finite( Pr_POP) & (Pr_POP >= 0) & (Pr_POP <= 1))) {
return( -1e100) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  ## The result!
  lglk <- sum( dpois( n_POP, lambda=n_comp_POP * Pr_POP, log=TRUE))

return( lglk)
}



"generic_lglk_fish" <-
function( params, log_prior=NULL, cheat=NULL){
  ## For debugging only:
  last_params <<- params

  ## Unpack parameters
  Ncad_y0 <- exp( params[ 1])           # total CATCHABLE abundance in year 1
  aslope_y0 <- exp( params[ 2])
  M <- exp( params[ 3])
  wtfecpar <- exp( params[ 4])
  fecata <- (weight ^ wtfecpar) * (1:AMAX >= Amat)


  ## NB constant recruitment !
  ## A "real" CKMR-commercial-fish model would also have Recruitment Deviations, one per cohort
  # ... implemented as Random Effects, with variance estimated PROPERLY eg using TMB / Laplace Approx
  # ... or MCMC, I suppose, if you must...
  # You would NOT use R
  # Quite a few other simplifications here compared to normal stock assessment...
  # ... mainly for clarity-of-exposition, also speed

  ## Population dynamics. No plus group --- certain death at AMAX
  # Just FEMALES
  # NB starts from first age where catches are taken--- no point in modelling younger ages explicitly
  N_ya <- offarray( 0, first=c( min_Bju, AminC), last=c( max_Bju, AMAX))

  # Split up numbers-at-age in first year
  init <- exp( -aslope_y0 * ALL_AGES) # actually it's not *all* ages, but AminC up
  N_ya[ y0,] <- Ncad_y0 * init / sum( init)
  N_ya[ ,AminC] <- c( N_ya[ y0, AminC]) # constant recruitment for all subsequent years

  if( any( is.na( N_ya))) {
return( SIGNAL_BADFIT) # MASSIVE negative number (setup in boring_data_prep):
    # ... Inf can cause problems, due to flaw in optim/nlminb
  }

  # Fill in all the pop dyn numbers
  # Should use Baranov (and NB M is pretty high for notog1), but Pope's equation is easiest
  eM2 <- exp( -M/2)
  for( y in tail( Bju_range, -1)) {
    N_ya[ y, Anonplus+1] <- (N_ya[ y-1, Anonplus]*eM2 - Cfem_ya[ y-1, Anonplus]) * eM2
  }

  if( any( N_ya < 0)) {
return( SIGNAL_BADFIT)
  }

  if( !is.null( cheat)){ # Normally this does nothing! But can use to
    # ... overwrite N_ya and fecata--- only for debugging!
    extract.named( cheat)
    N_ya <- N_sya[ SLICE='F', Bju_range, ALL_AGES] # keep the bits we want
    fecata <- fec_sa[ SLICE='F', 1:AMAX]
  }


  ## CKMR stuff goes here...
  TRO <- offarray(
      N_ya[,ADAGES] %**% fecata[ ADAGES],
      first=min_Bju,
      last=max_Bju)

  # Define this here so it knows about fecata. Vectorized (ie multiple 'a' at once).
  fec_fn <- function( a) {
      # check 0 <= a < AMAX
      # If not, change a to something that allows fecata[ awork] to succeed...
      awork <- pmax( 1, pmin( AMAX, a))
      # ... and set fec to 0 if a <= 0
      fecata[ awork] * (a > 0)
    }


  Pr_POP <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=ADAGES, {
      Bad <- Yad - Aad
      ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
      ( Bju >= Amat + Bad ) *   # was ad mature at B[ju] ?
      ( fec_fn( Aad - (Yad-Bju)) / TRO[ Bju]) # ERRO if above
    })

  # Cumul surv
  psurv_ayy <- autoloop( a1=ADAGES, y1=Bju_range, y2=Bju_range, {
      dy <- y2-y1
      a2 <- a1 + abs( dy) # dy<0 will be zapped below--- this is anti-overrun
      a2_clip <- pmin( a2, AMAX) # avoid index overrun

      ans <- (dy >= 0) *      # forwards only
        (a2 <= AMAX) *   # last age class dies entirely
        ( N_ya[ y2, a2_clip] / N_ya[ y1, a1]) # THIS is most important line
    return( ans)
    })

  # Cond age of a Mother at B1: do-able in basic R, without autoloopery etc
  Pr_Mum1Age <- autoloop( y=Bju_range, a=ADAGES, {
    return( N_ya[ y, a] * fecata[ a] / TRO[ y])
    })

  # Prob of HSP with Mum also being of specific age when first was born
  Pr_HSP_and_Mum1Age <- autoloop( B1=Bju_range, B2=Bju_range, PARAGE=ADAGES, {
      # Conditional on Mum's age at B1 (ie PARAGE)
      Pr_HSP_cond <- (B1 < B2) * psurv_ayy[ PARAGE, B1, B2] * fec_fn( PARAGE+B2-B1) / TRO[ B2]

      # NB indices [B1,B2,PARAGE] are implicit in next line
    return( Pr_HSP_cond * Pr_Mum1Age[ B1, PARAGE])
    })

  # That's the joint prob of HSPness and Mum-age. We just want the marginal prob of HSPness
  Pr_HSP <- sumover( Pr_HSP_and_Mum1Age, 'PARAGE')
  Pr_HSP <- Pr_HSP * (1-Pr_FNeg_HSP) # allow for known loss-rate due to kinference

  lglk <- sum( dpois( n_POP, size=n_comp_POP, prob=Pr_POP, log=TRUE)) +
      sum( dpois( n_HSP, size=n_comp_HSP, prob=Pr_HSP, log=TRUE))

  # Optional penalty on the params, to avoid getting stuck...

  if( !is.null( log_prior)) {
    lglk <- lglk + log_prior( params)
  }


  if( !is.null( env$verbose)) { # for debugging
    print( params)
    scatn( 'HSP O & E: %i & %5.2f', sum( n_HSP), sum( n_comp_HSP * Pr_HSP))
    scatn( 'POP O & E: %i & %5.2f', sum( n_POP), sum( n_comp_POP * Pr_POP))
    if( sum( n_HSP_sc) > 0) { # check same-cht HSPs, if they are collected; in this sim, they aren't...
      # ... and they won't (shouldn't) be used directly in estimation, either
      scatn( 'HSP same-cohort O & E FYI:  %i & %5.2f',
          sum( n_HSP_sc), sum( n_comp_HSP_sc * diag( unclass( Pr_HSP))))
    }
  }


  # Useful to keep some stuff after function exits--- you can
  # ... get this from "env" afterwards. Trust me, this works...
  list2env( mget( cq( Ncad_y0, M, fecata, TRO, N_ya, Pr_POP, Pr_HSP)), env)

return( lglk)
}



"generic_lglk_HSP_mammal" <-
function( params){
  ## Unpack parameters
  Nfad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  Z <- exp( params[ 3])

  ## Population dynamics

  N <- offarray( 0, first=min( Bju_range), last=max( Bju_range))
  N[] <- Nfad_y0 * exp( RoI * (Bju_range-y0))
  # or:
  # N <- autoloop( B=Bju_range, {
  #   Nfad_y0 * exp( RoI * (B-y0)) 
  # })

  # OMHSP is "Ordered" Maternal HSP--- ie 1st one was born first
  # Expects data is organized that way (and thus that we know age exactly)
  # Watch out for double-counting with HSPs...
  Pr_OMHSP <- autoloop( B1=Bju_range, B2=Bju_range, {
    cumul_surv <- exp( -Z*(B2-B1))
    (B2 > B1) * cumul_surv / N[ B2]
  })

  Pr_MHSP <- Pr_OMHSP # for consistency
 
  # Useful to keep some stuff after function exits. Trust me, this works...
  list2env( mget( cq( Nfad_y0, RoI, Z, N, Pr_MHSP)), env)

  if( !all( is.finite( Pr_MHSP) & (Pr_MHSP >= 0) & (Pr_MHSP <= 1))) {
return( -1e100) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  lglk <- sum( dpois( n_MHSP, lambda= n_comp_MHSP * Pr_MHSP, log=TRUE))

return( lglk)
}



"generic_lglk_MOP_fuzzage_adult_mammal" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  # ... that's all,  folks

  ## Population dynamics
  N <- offarray( 0, first=min( years), last=max( years))
  N[] <- Nfad_y0 * exp( RoI * (years-y0))

  ## Ideal kinprob

  Pr_MOP_ideal_BYA <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=likely_Aad_range, {
    Bad <- Yad - Aad
    ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
    ( Bju >= Amat + Bad ) *   # was ad mature at B[ju] ?
    ( 1/N[ Bju])                  # ad's ERO at B[ju] if alive & mature, divided by TRO that year
  })

  ## Kinprob given OBSERVED covars, ie Fuzzage instead of true A
  
  ## Bayes' theorem to get Pr[ truA | Fuzzage, Y]
  Pr_Afuzz_Y <- autoloop( A=likely_Aad_range, Fuzz=Fuzzad_range, Y=Yad_range, {
    Pr_Fuzzage_Atru[ Fuzz, A] * 
    est_Pr_A_SampY[ A, Y]
  })
  
  Pr_Fuzz_Y <- sumover( Pr_Afuzz_Y, 'A')
  Pr_A_FuzzY <- autoloop( A=likely_Aad_range, Fuzz=Fuzzad_range, Y=Yad_range, {
    Pr_Afuzz_Y[ A, Fuzz, Y] / Pr_Fuzz_Y[ Fuzz, Y]
  })
    
  # "Mix" the ideal probs 
  Pr_AMOP_BYFuzz <- autoloop( Bju=Bju_range, Y=Yad_range, A=likely_Aad_range, Fuzz=Fuzzad_range, {
    Pr_MOP_ideal_BYA[ Bju, Y, A] *
    Pr_A_FuzzY[ A, Fuzz, Y]
  })
  Pr_MOP_BYFuzz<- sumover( Pr_AMOP_BYFuzz, 'A')

  # Check: dimseq( Pr_MOP_BYFuzz) vs dimseq( n_MOP)
  Pr_MOP <- Pr_MOP_BYFuzz #  keep same notation as other examples

  # Useful to keep some stuff after function exits. Trust me, this works...
  list2env( mget( cq( Nfad_y0, RoI, N, Pr_MOP)), envir=environment( sys.function()))

  if( !all( is.finite( Pr_MOP) & (Pr_MOP >= 0) & (Pr_MOP <= 1))) {
return( -1e100) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  lglk <- sum( dpois( n_MOP, lambda=n_comp_MOP * Pr_MOP, log=TRUE))

return( lglk)
}



"generic_lglk_MOP_fuzzall_mammal" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  # ... that's all,  folks

  ## Population dynamics
  N <- autoloop( Y=years, {
    Nfad_y0 * exp( RoI * (Y-y0))
  })

  ## Ideal kinprob
  Pr_MOP_BYA <- autoloop( Bju=years, Ypar=Y_range, Apar=A_range, {
    Bpar <- Ypar - Apar
    ( Bju <= Ypar ) *          # was ad still alive at B[ju] ?
    ( Bju >= Amat + Bpar ) *   # was ad mature at B[ju] ?
    ( 1/N[ Bju])                  # ad's ERO at B[ju] if alive & mature, divided by TRO that year
  })

  ## Kinprob given OBSERVED covars, ie Fuzzage instead of true A
  
  ## Bayes' theorem to get Pr[ truA | Fuzzage, Y]
  Pr_AF_Y <- autoloop( A=A_range, F=F_range, Y=Y_range, {
    Pr_Fuzzage_Atru[ F, A] * 
    est_Pr_A_SampY[ A, Y]
  })
  
  Pr_F_Y <- sumover( Pr_AF_Y, 'A')
  Pr_A_FY <- autoloop( A=A_range, F=F_range, Y=Y_range, {
    Pr_AF_Y[ A, F, Y] / Pr_F_Y[ F, Y]
  })
    
  # "Mix" the ideal probs 
	Pr_MOPAA_FYFY <- autoloop( A1=A_range, A2=A_range, F1=F_range, Y1=Y_range, F2=F_range, Y2=Y_range, {
		# True kinship (MO or OM) is unobservable. MO and OM are mutually exclusive, so we can add their probs...
		# NB horrible truncation here... have to avoid accessing OOR elements. 
		# If one animal is clearly adult, the truncation is "bad" but the prob of it being an off anyway will be 0
		( Pr_MOP_BYA[ pmax( y0, Y1-A1), Y2, A2] + Pr_MOP_BYA[ pmax( y0, Y2-A2), Y1, A1]) *    
		 Pr_A_FY[ A1, F1, Y1] *
		 Pr_A_FY[ A2, F2, Y2]
	})

	Pr_MOP_FYFY <- sumover( Pr_MOPAA_FYFY, cq( A1, A2))

	rm( Pr_MOPAA_FYFY); gc() # _might_ speed things up...

  # Check: dimseq( Pr_MOP_FYFY) vs dimseq( n_MOP)
  Pr_MOP <- Pr_MOP_FYFY #  keep same notation as other examples

  # Useful to keep some stuff after function exits. Trust me, this works...
  list2env( mget( cq( Nfad_y0, RoI, N, Pr_MOP)), envir=environment( sys.function()))

  if( !all( is.finite( Pr_MOP) & (Pr_MOP >= 0) & (Pr_MOP <= 1))) {
return( -1e100) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  lglk <- sum( dpois( n_MOP, lambda=n_comp_MOP * Pr_MOP, log=TRUE))

return( lglk)
}



"generic_lglk_MOP_ideal_mammal" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  # ... that's all,  folks

  ## Population dynamics

  N <- offarray( 0, first=min( years), last=max( years))
  N[] <- Nfad_y0 * exp( RoI * (years-y0))

  # offarray::autoloop() evaluates an expression over all combos of the args
  # ... and returns an array indexed by all combos. Like using nested loops.
  # Indices must match n_comp[] and n_MOP[], which are set up by boring_data_prep()

  Pr_MOP <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=Aad_range, {
    Bad <- Yad - Aad
    ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
    ( Bju >= Amat + Bad ) *   # was ad mature at B[ju] ?
    ( 1/N[ Bju])                  # ad's ERO at B[ju] if alive & mature, divided by TRO that year
  })

  # Useful to keep some stuff after function exits. Trust me, this works...
  list2env( mget( cq( Nfad_y0, RoI, N, Pr_MOP)), envir=environment( sys.function()))

  if( !all( is.finite( Pr_MOP) & (Pr_MOP >= 0) & (Pr_MOP <= 1))) {
return( -1e100) # optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dbinom()
  }

  lglk <- sum( dpois( n_MOP, lambda=n_comp_MOP * Pr_MOP, log=TRUE))

return( lglk)
}



"generic_lglk_MOP_LR_mammal" <-
function( params){
  assign( 'last_params', params, environment( sys.function())) # debugging etc

  ## Unpack parameters
  Nfad_y0 <- exp( params[ 1])
  RoI <- params[ 2]
  p <- inv.logit( params[3]) # prob of r-gene; R-pheno <=> 2 r-genes

  ## Population dynamics

  N <- offarray( 0, first=min( years), last=max( years))
  N[] <- Nfad_y0 * exp( RoI * (years-y0))

  # offarray::autoloop() evaluates an expression over all combos of the args
  # ... and returns an array indexed by all combos. Like using nested loops.
  # Indices must match n_comp[] and n_MOP[], which are set up by boring_data_prep()

  HANDED <- cq( L, R)

  Pr_H <- c( R=sqr( p), L=1-sqr(p))
  
  Pr_Hj_HaMO <- offarray( 0, first=c(1,1), last=c(2,2), dimnames=list( Ha=HANDED, Hj=HANDED))
  Pr_Hj_HaMO[ 'L', 'L'] <- 1 - sqr(p)/(1+p)
  Pr_Hj_HaMO[ 'R', 'L'] <- sqr(p)/(1+p)
  Pr_Hj_HaMO[ 'L', 'R'] <- 1-p
  Pr_Hj_HaMO[ 'R', 'R'] <- p
  
  # Conditional ERRO (ie Pr[MO]), given ad was alive and mature, and the handednesses. See maths doco
  cond_ERRO_almatBHH <- autoloop( Bju=Bju_range, Hju=HANDED, Had=HANDED, {
    (Pr_Hj_HaMO[ Hju, Had] / Pr_H[ Hju]) * (1/N[ Bju])               
  })
  
  Pr_MOP <- autoloop( Bju=Bju_range, Hju=HANDED, Yad=Yad_range, Aad=Aad_range, Had=HANDED, {
    Bad <- Yad - Aad          # when was ad born?
    ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
    ( Bju >= Amat + Bad ) *   # was ad mature at B[ju] ?
    cond_ERRO_almatBHH[ Bju, Hju, Had] 
  })

  # Useful to keep some stuff after function exits. Trust me, this works...
  list2env( mget( cq( Nfad_y0, RoI, p, N, Pr_H, Pr_MOP)), envir=environment( sys.function()))

  if( !all( is.finite( Pr_MOP) & (Pr_MOP >= 0) & (Pr_MOP <= 1))) {
return( -1e100) # ie optim/nlminb has tried insane param values...
    # ... this avoids NaN warnings from dpois()
  }

  lglk <- sum( dpois( n_MOP, lambda= n_comp_MOP * Pr_MOP, log=TRUE))

return( lglk)
}



"M_fun" <-
function( par) { lglk_fish( par); with( env,
    M
  )}



"pairid" <-
function( i, j) i + 1/(1+j)



"pyro_enfum" <-
function( AMIN, AMAX, mulbo=1) {
  ## Bigger mulbo means tighter
  truA <- AMIN:AMAX
  
  cumpr_aconv_a <- autoloop( Aconv=0:AMAX, Atru=truA, 
    pgamma( Aconv, shape=Atru*mulbo, scale=1/mulbo)
  )

  pr_aconv_a <- autoloop( Aconv=1:AMAX, Atru=truA, 
    cumpr_aconv_a[ Aconv, Atru] - cumpr_aconv_a[ Aconv-1, Atru]
  )
    
  # Gotta get SOME recorded age for each true age, so norm this...
  flurk <- sumover( pr_aconv_a, 'Aconv')
  pr_aconv_a <- autoloop( Aconv=1:AMAX, Atru=truA, 
    pr_aconv_a[ Aconv, Atru] / flurk[ Atru]
  )

return( pr_aconv_a)
}



"sample_rhombi" <-
function( m_ju, m_ad, sel_m, Nm, Nf, fec_Co, sel_Co, Pr_Co) {
	Co_levels <- seq_along( fec_Co)
	n_Co_levels <- length( fec_Co)
stopifnot( length( sel_Co)==n_Co_levels)	

  Co_m <- rsample( Nm, Co_levels, prob=Pr_Co, replace=TRUE)

  # Only need ancestry for the samples
  Dad <- rsample( m_ju, seq_len( Nm), prob=fec_Co[ Co_m], replace=TRUE)
  Mum <- rsample( m_ju, seq_len( Nf), replace=TRUE)
  
  # Derwent: unselective sampling
  msamp_m <- floor( m_ad * (sel_m * Nm ) / (sel_m * Nm + Nf))
  ad_m <- rsample( msamp_m, seq_len( Nm), replace=FALSE, prob=sel_Co[ Co_m])
  msamp_f <- m_ad - msamp_m
  ad_f <- rsample( msamp_f, seq_len( Nf), replace=FALSE)
  
  # Women & children first! So, offset the ids of Juves by msamp_f, and of Male adults by msamp_f + m_ju
  iMOPs <- match( Mum, ad_f, 0)
  MOPs <- cbind( iMOPs %except% 0, msamp_f + which( iMOPs > 0))
  
  iFOPs <- match( Dad, ad_m, 0)
	FOPs <- cbind( msamp_f + m_ju + (iFOPs %except% 0), msamp_f + which( iFOPs > 0))
	POPs <- rbind( MOPs, FOPs)
	POPs <- POPs[ rsample( nrow( POPs), 1:nrow( POPs), replace=TRUE),]

  m <- m_ad + m_ju
  LETTAZ <- matrix( rsample( 5*m, LETTERS, replace=T), ncol=5)
  random_names <- LETTAZ[,1] %&% LETTAZ[,2] %&% LETTAZ[,3] %&% LETTAZ[,4] %&% LETTAZ[,5]
  
  Samps <- data.frame(
    Me=random_names,
    Sex=c( rep( 'F', msamp_f), rep( NA, m_ju), rep( 'M', msamp_m)),
    Co=c( rep( 1, msamp_f), rep( 0, m_ju), Co_m[ ad_m]),
    poss_par= c( rep( TRUE, msamp_f), rep( FALSE, m_ju), rep( TRUE, msamp_m)),
    poss_off= c( rep( FALSE, msamp_f), rep( TRUE, m_ju), rep( FALSE, msamp_m))
  )

  stuff <- returnList( Samps, POPs)
  stuff@public <- mget( cq( Co_levels, n_Co_levels))
return( stuff)
}



"SPRR_fun" <-
function( par, ranges=FALSE) {
  lglk_fish( par)
  env$ranges <<- ranges
  result <- with( env, {
    # Reconstruct F...
    lastNyears <- head( tail( dimseq( N_ya)[[1]], 3), -1)
    relevA <- 4 %upto% (dimrange( N_ya)[2,2]-1)
    if( ranges) { # can't test directly...
  return( returnList( lastNyears, relevA))
    }

    Zlast <- -log( N_ya[ lastNyears+1, relevA+1] / N_ya[ lastNyears, relevA])
    Flast <- Zlast - M
    Fav <- colMeans( Flast)
    Zav <- colMeans( Zlast)
    Mav <- 0*Zav + M # KISS!
    Na_now <- exp( -cumsum( c( 0, Zav)))
    Na_virgin <- exp( -cumsum( c( 0, Mav)))
    relevFec <- fecata[ 4:(dimrange( N_ya)[2,2])]
    ROlife_virgin <- Na_virgin %*% relevFec
    ROlife_now <- Na_now %*% relevFec
    SPRR <- ROlife_now / ROlife_virgin
  return( c( SPRR))
  })
  env$ranges <<- NULL
return( result)
}



"TRO7_est_fun" <-
function( pars) { lglk_fish( pars); with( env,
    N_ya %**% fecata[ AFEC] / fecata[SLICE=7]
  )}



"vbdiff_fec" <-
function( k, n=100, Amax=20, Amat=3, percent_sd=15, dby=0.2) {
  qq <- (1:n - 0.5) / n
  Linf <- qnorm( qq, 100, sd=percent_sd)

  lifespan <- pmin( pmax( rexp( n, 3/Amax), 1), Amax)

  A <- 1:Amax

  require( vecless) # easier than outer() at least for me...


  vbl <- function( linf, a) linf * (1-exp( -k*a))
  L[ l, a]:= vbl( Linf[ l], A[a])
  # L <- outer( Linf, A, vbl) # outer is confusing

  plot( 0, 0, xlim=range( A), ylim=c( 0, max( Linf)), type='n',
      xlab='Age', ylab='Length', cex.lab=2, cex.axis=2, mar=c( 5, 5, 2, 1)+0.1)
  abline( v=Amat, col='grey')
  text( Amat, 100, 'Maturity', adj=c( 1, 0), cex=1.5, srt=90, col='grey')
  # matplot is shite
  for( i in 1:n) {
    afrac <- seq( 1, lifespan[ i], by=dby)
    lines( afrac, vbl( Linf[ i], afrac))
  }

  fec[ l, a] := L[ l, a]^3 * (lifespan[ l] >= A[ a])
  # Fec vblty

  Pr_HSP_true <- sum( fec*fec) / sqr( sum( fec))

  fecabar[ a] := (SUM_ %[l]% (L[ l, a]^3) ) / n

  nata[ a]:= SUM_ %[l]% (lifespan[ l] >= A[ a])
  Pr_HSP_justa <- (nata %**% sqr( fecabar)) / sqr( nata %**% fecabar)

return( c( bias_Nhat = Pr_HSP_justa / Pr_HSP_true))
}




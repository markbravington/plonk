function( params){


  # Remember to NOT EVEN DO comparisons when Hank caught in spawning season of Irene's birth



  for( B in all_birthyears)
    for( S in {MALES, FEMALES})
      # work out TRO
      temp <- 0;
      for( age=1 to maxage) { 
        temp <- temp + N[ B, S, A] * fec[ A, S]
      }
      TRO[ S, B ] <- temp
    }
  }


  for( all HANKs)
  for( all_IRENEs) {
   Pr_K_HANK_IRENE_FO <- 
     if( b[ IRENE] >= y[ HANK]) = 0 else 
     if( b[ HANK] >= b[ IRENE] 0 else
       fec[ Male, Hanks age at Irene's birth] / TRO[ Male, Irene's birth] 
     };
   
   if( K_is_FO[ HANK, IRENE]) {
     lglk <- lglk + log( Pr_K_HANK_IRENE_FO) 
   } else {
     lglk <- lglk + log( 1-Pr_K_HANK_IRENE_FO)
   }
  }

# Version from lglk_fish	

  fec_fn <- function( a) {
      # check 0 <= a < AMAX
      # If not, change a to something that allows fecata[ awork] to succeed...
      awork <- pmax( 1, pmin( AMAX, a))
      # ... and set fec to 0 if a <= 0
      fecata[ awork] * (a > 0)
    }


  Pr_POP <- autoloop( Bju=Bju_range, Yad=Yad_range, Aad=ADAGES, {
      Bad <- Yad - Aad;
      ( Bju <= Yad ) *          # was ad still alive at B[ju] ?
      ( Bju >= Amat + Bad ) *   # was ad mature at B[ju] ?
      ( fec_fn( Aad - (Yad-Bju)) / TRO[ Bju]) # ERRO if above
    })


  for( SAM...)
    for( HOLLY...) {
      if( b[ HOLLY] < b[ SAM]){
        Pr_K_SAM_HOLLY_is_OMHSP <- 0
      }
      temp <- 0
      for( A = 1 to 20) {
        temp <- temp + (fec_fn( A, FEMALE)/TRO[] * Pr_surv[ A, b[SAM], b[HOLLY]] 8 
      }
      Pr_K_SAM_HOLLY_is_OMHSP <- temp

      if( 


}

fe
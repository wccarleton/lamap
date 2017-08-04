#' LAMAP main functions
#'
#' This function
#'
#'
#'
#'
lamap <- function(observed,
                  knownsite_pcdfs,
                  knownsite_coords,
                  steps,
                  maxsites=NULL,
                  weightfun=NULL,
                  weightparams=NULL,
                  combinations=NULL,
                  ...){
   if(ncol(observed) < 3){
      stop(paste("Missing columns. There should be two leading columns ",
                 "in observed for x and y coordinates, followed by one ",
                 "or more variable columns.",
                 sep=""))
   }
   sitedists <- orderSites(test_observed,knownsite_coords)
   if(!is.null(maxsites)){
      sites_included <- sitedists[1:maxsites,"index"]
   } else {
      sites_included <- sitedists[,"index"]
   }
   n_knownsites <- length(knownsite_pcdfs[sites_included])
   if(is.null(combinations)){
      combinations <- prepCombinations(n_knownsites)
   }
   jdensities <- lapply(knownsite_pcdfs[sites_included],function(jds){
      jdensity(observed=observed[,-c(1,2)],steps=steps,pcdfs=jds,...)
   })
   prob_scores <- sapply(jdensities,function(x){
      prod(x[,2]-x[,1])
   })
   if(!is.null(weightfun)){
      weights <- weight(sitedists[1:maxsites,"dist"],weightfun,weightparams)
      prob_scores = prob_scores * weights
   }
   prob_union <- unionIndependent(prob_scores,combinations=NULL)
   return(prob_union)
}

#' lamap
#'
#' The primary lamap function. It estimates a lamap value given a single
#' observation and a training dataset based on known archaeological sites. The
#' observation is a single raster cell that might have multiple variable
#' dimensions—i.e., a single raster cell from a raster stack. Most analyses
#' will need to be parallelized. This is the function that should be called in
#' in parallel. It can be thought of as a computational kernel that processes a
#' single raster cell from the desired output raster map (i.e., the lamap surface).
#' The lamap package contains a wrapper function for parallel processing using
#' the R package 'snow'—see 'parLamap'.
#'
#' @param observed Single row dataframe of observed data (will be coerced).
#'  Each column should be one variable that corresponds to one of the cdfs in
#'  pcdfobj.
#' @param knownsite_pcdfs List of lists containing cdfs for each known site.
#'  See 'pcdf()'. This parameter should be the output of the function
#'  'knownsitePcdfs'.
#' @param knownsite_coords Dataframe containing xy coordinates for each known
#'  site. The columns need to be named "id","x","y". The ids should be ordered
#'  to correspond to the order of sites in the knownsite_pcdfs list. For now,
#'  only UTM coordinates should be used. See 'knownsiteCoords', a convenience
#'  function for creating this list from a dataframe.
#' @param steps A vector of window sizes for each variable of interest. The
#'  steps will be used to estimate integrals from the knownsite_pcdfs. They
#'  will be added/subtracted from the observed values to create integration
#'  intervals. The order of steps in this vector must correspond to the order
#'  of column variables in the 'observed' dataframe—which must also match the
#'  order of columns in training dataframe—i.e., 'traindf'—passed to the 'pcdf'
#'  function.
#' @param maxsites The maximum number of sites from the known-site training data
#'  to include in the calculations. This is primarily a convenience parameter
#'  that can be used to speed up calculations when the training set contains
#'  a large number of sites.
#' @param weightfun The function used to weight the lamap values by distance.
#'  It can be one of 'uniform' or 'exponential'—see 'weight' function. Default
#'  is Null, meaning no weighting is applied.
#' @param weightparams A vector of parameters passed to the 'weight' function.
#'  Default is Null.
#' @param combinations A conveience paramter used to speed up calculations. The
#'  lamap calculation involves the use of hte inclusion-exclusion principle for
#'  estimating the total probability of a union of independent events.
#'  Consequently, the R 'combn' function is used and can slow the calculations
#'  considerably. However, the user can instead use a utility function provided
#'  with this package—'prepCombinations'. It computes the combinations which can
#'  then be saved and passed to the lamap function with this parameter, saving
#'  the algorithm from having to regenerate the combinations for every cell in
#'  in the desired output map. Default is Null.
#' @param ... Additional paramters passed to the 'jdensity' function:
#'  the 'nosupport', 'partial', and 'interpolate' options.
#' @return A single lamap value.
#' @export

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
   sitedists <- orderSites(observed,knownsite_coords)
   if(!is.null(maxsites)){
      sites_included <- sitedists[1:maxsites,"index"]
   } else {
      sites_included <- sitedists[,"index"]
   }
   n_knownsites <- length(knownsite_pcdfs[sites_included])
   if(is.null(combinations)){
      combinations <- prepCombinations(n_knownsites)
   }
   if(any(is.na(observed))){
      prob_union <- NA
   } else {
      jdensities <- lapply(knownsite_pcdfs[sites_included],function(jds){
         jdensity(observed=observed[,-c(1,2)],steps=steps,pcdfs=jds,...)
         })
      prob_scores <- sapply(jdensities,function(x){
         prod(x[,2]-x[,1])
         })
      if(!is.null(weightfun)){
         weights <- weight(sitedists[1:length(prob_scores),"dist"],weightfun,weightparams)
         prob_scores = prob_scores * weights
      }
      prob_union <- unionIndependent(prob_scores,combinations)
   }
   return(prob_union)
}

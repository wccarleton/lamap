#' knownsiteCoords
#'
#' Calculates the probability of a union of independent events using the
#' Law of Total probability and the inclusion-exclusion principle.
#'
#' @param knownsitedf A dataframe containing raster data from known sites.
#' @return A dataframe of sites with centre coordinates.
#' @export

knownsiteCoords <- function(knownsitedf){
   siteid <- unique(knownsitedf$id)
   if(is.null(siteid)){
      stop("Error: missing site ids")
   }
   if(is.null(knownsitedf$x) | is.null(knownsitedf$y)){
      stop("Error: missing site coordinates")
   }
   site_coords <- lapply(siteid,function(j,knownsitedf){
      site_data <- knownsitedf[which(knownsitedf$id==j),]
      sitecoord_x = mean(range(site_data$x))
      sitecoord_y = mean(range(site_data$y))
      cbind(sitecoord_x,sitecoord_y)
   },knownsitedf)
   site_coords <- do.call(rbind,site_coords)
   knownsite_coords <- cbind(siteid,site_coords)
   return(knownsite_coords)
}

#' knownsitePcdfs
#'
#' Calculates the probability of a union of independent events using the
#' Law of Total probability and the inclusion-exclusion principle.
#'
#' @param probs A vector of probabilities
#' @param combinations A list of matrices containing the the possible
#'  of probabilities. This values in the matrix cells should be the indeces of
#'  the probabilities in the vector 'probs'. Default is Null—so the combinations
#'  will be arranged by the function. To create this combination list, you can
#'  use the 'prepCombinations' convenience function provided in the 'lamap'
#'  package.
#' @return Probability of the union.
#' @export

knownsitePcdfs <- function(knownsitedf){
   if(ncol(knownsitedf) <= 3){
      stop("Missing columns. There should be three leading columns for id, x, and y, followed by one or more variable columns.")
   }
   siteid <- unique(knownsitedf$id)
   if(is.null(siteid)){
      stop("Error: missing site ids")
   }
   knownsite_pcdfs <- lapply(unique(knownsitedf$id),function(x){
      pcdf(knownsitedf[which(knownsitedf$id==x),-c(1:3)],interpolate=F)
      })
   return(knownsite_pcdfs)
}

#' numericFactorLevels
#'
#' Returns numeric factor levels from a Table Factor—used internally.
#'
#' @param x Table Factor.
#' @return Numeric factor levels.
#' @export

numericFactorLevels <- function(x){
   return(cbind(as.numeric(levels(x[,1])),x[,2]))
}

#' findOverlap
#'
#' Determines the relationship between the bounds of two sets (ranges). Imagine
#' two sets, 'a' and 'b', where 'a' has bounds represented graphically by []
#' and 'b' has bounds ()—note this is not intended to represent set notation.
#' The function returns 1 if 'b' is entirely within 'a', like this: [()]. It
#' return 2 if [(]); 3 if [](); 4 if ([]); 5 if ([)]; and lastly 6 if ()[]. This
#' function is used internally if the 'overlap' paramater is True in the 'lamap'
#' function.
#'
#' @param a A vector of 2 dimensions—e.g., (1,10)
#' @param b A vector of 2 dimensions—e.g., (5,9)
#' @return Int. Indicating relationship of boundaries.
#' @export

findOverlap <- function(a,b){
   #[a](b)
   if(a[1] < b[1] & a[2] > b[2]){
      #[()]
      return(1)
   } else if(a[1] < b[1] & a[2] <= b[2] & a[2] >= b[1]){
      #[(])
      return(2)
   } else if(a[2] < b[1]){
      #[]()
      return(3)
   } else if(a[1] >= b[1] & a[2] <= b[2]){
      #([])
      return(4)
   } else if(a[1] >= b[1] & a[1] < b[2] & a[2] >= b[2]){
      #([)]
      return(5)
   } else if(a[1] >= b[2]){
      #()[]
      return(6)
   }
}

#' orderSites
#'
#' Orders known sites (from the knownsite_coords dataframe) according to their
#' straight-line (Euclidean) distance from a given observed raster cell. It is
#' used when calculating and applying weights to the LAMAP values.
#'
#' @param observed A single row dataframe with at least one column named 'x' and
#'  one named 'y', corresponding to the xy coordinates of a raster cell in the
#'  'rasterdata' raster object.
#' @param knownsite_coords A dataframe constucted with the 'knownsiteCoords'
#'  function, or a dataframe with at least three columns named 'id','x', and 'y'
#' @return A dataframe with a column of sorted distances.
#' @export

orderSites <- function(observed,knownsite_coords){
   observed_xy <- observed[,c("x","y")]
   knownsite_xy <- knownsite_coords[,c("x","y")]
   dists <- pointDistance(observed_xy,knownsite_xy,lonlat=F)
   knownsite_dists <- data.frame(id=knownsite_coords$id,dist=dists)
   theorder <- order(knownsite_dists[,2])
   knownsite_dists_sorted <- cbind(knownsite_dists[theorder,],index=theorder)
   return(knownsite_dists_sorted)
}

#' prepCombinations
#'
#' Prepares a list of combinations (matrices containing the possible
#' combinations of n).
#'
#' @param n The number of things to combine in 1–n ways.
#' @return A length n list of matrices.
#' @export

prepCombinations <- function(n){
   combinations <- sapply(1:n,function(x,l)combn(c(1:l),x),n)
   return(combinations)
}

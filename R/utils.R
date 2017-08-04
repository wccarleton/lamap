#' LAMAP utility functions
#'
#' These functions provide support for manipulating data.
#'
#'
#'
#'
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

knownsitePcdfs <- function(knownsitedf){
   if(ncol(knownsitedf) <= 3){
      stop("Missing columns. There should be three leading columns for id, x, and y, followed by one or more variable columns.")
   }
   siteid <- unique(knownsitedf$id)
   if(is.null(siteid)){
      stop("Error: missing site ids")
   }
   knownsite_pcdfs <- lapply(unique(xyz_nona_rnd$id),function(x){
      pcdf(xyz_nona_rnd[which(xyz_nona_rnd$id==x),-c(1:3)],interpolate=F)
      })
   return(knownsite_pcdfs)
}

numericFactorLevels <- function(x){
   return(cbind(as.numeric(levels(x[,1])),x[,2]))
}

findClosest <- function(x,v){
   insupport <- (x >= min(v) && x <= max(v))
   if(insupport){
      if(any(x == v)){
         return(x)
      } else {
         lower_index <- which(diff(v > x) == 1)
         upper_index = lower_index + 1
         lower_diff = abs(x - v[lower_index])
         upper_diff = abs(x - v[upper_index])
         if(lower_diff <= upper_diff){
            return(v[lower_index])
            } else {
               return(v[upper_index])
               }
         }
   } else {
      return(NA)
   }
}

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

utmDist <- function(x1,x2){
   dist = sqrt((x1[1]-x2[1])^2+(x1[2]-x2[2])^2)
   return(dist)
}

orderSites <- function(observed,knownsite_coords){
   observed_xy <- observed[,c("x","y")]
   knownsite_xy <- knownsite_coords[,c("x","y")]
   dists <- pointDistance(observed_xy,knownsite_xy,lonlat=F)
   knownsite_dists <- data.frame(id=knownsite_coords$id,dist=dists)
   theorder <- order(knownsite_dists[,2])
   knownsite_dists_sorted <- cbind(knownsite_dists[theorder,],index=theorder)
   return(knownsite_dists_sorted)
}

prepCombinations <- function(n){
   combinations <- sapply(1:n,function(x,l)combn(c(1:l),x),n)
   return(combinations)
}

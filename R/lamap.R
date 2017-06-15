#' LAMAP main function
#'
#' This function
#'
#'
#'
#'
lamap <- function(observed,knownsites,...){

}

knownsite <- function(knownsitedf){
   nvars <- ncols(knownsitedf)
   siteid <- unique(knownsitedf$id)
   if(is.null(siteid) | length(siteid) > 1){
      stop("Error: missing or non-unique site id")
   }
   if(is.null(knownsitedf$x) | is.null(knownsitedf$y)){
      stop("Error: missing site coordinates")
   }
   sitecoord_x <- (max(knownsitedf$x) - min(knownsitedf$x))/2
   sitecoord_y <- (max(knownsitedf$y) - min(knownsitedf$y))/2
   pcdfs <- pcdf(knownsitedf[,c(4:nvars)])
   site <- list(siteid=siteid,
      sitecoords=c(sitecoord_x,sitecoord_y),
      pcdfs=pcdfs)
   return(site)
}

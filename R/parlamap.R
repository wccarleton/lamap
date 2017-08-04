parLamapCaller <- function(x,
                           rasterdata,
                           knownsite_pcdfs,
                           knownsite_coords,
                           steps,
                           maxsites,
                           nosupport,
                           partial) {
   cat(paste("\r",x))
   observedcell <- getObservedData(x, rasterdata)
   lamapoutput <- lamap(observed=observedcell,
                        knownsite_pcdfs=knownsite_pcdfs,
                        knownsite_coords=knownsite_coords,
                        steps=steps,
                        maxsites=maxsites,
                        nosupport=NA,
                        partial=T)
   return(lamapoutput)
}

getObservedData <- function(x, rasterdata){
   observedcell <- as.data.frame(cbind(xyFromCell(rasterdata, x),
                                 round(extract(rasterdata, x))))
   return(observedcell)
}

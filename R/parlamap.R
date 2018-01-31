#' parLamap
#'
#' Parallel wrapper for the primary 'lamap' function. The wrapper is based on
#' 'snow' clusters. It assumes that the cluster environment is such that the
#' nodes share storage space—i.e., each has direct access to the same folder
#' specified in the 'rasterpath' argument. It also assumes that the nodes have
#' the same R environments, including all lamap functions and the
#' knownsite_coords' and 'knownsite_pcdfs' objects. If your setup is different,
#' you will have to create your own parallel scripts.
#'
#' @param cluster_object A 'snow' cluster returned by a function like
#' 'makeCluster'—see the 'snow' package.
#' @param rasterpath String. Path to geoTiff stack containing spatial variables.
#'  It will also provide the spatial meta data for the output lamap surface
#'  (i.e., it will define the prediction region and spatial resolution).
#' @param outputpath String. Path to output folder.
#' @param knownsite_pcdfs List of of lists containing cdfs for each known site.
#'  See 'pcdf()'.
#' @param knownsite_coords Dataframe containing xy coordinates for each known
#'  site. The columns need to be named "id","x","y". The ids must be ordered
#'  to correspond to the order of sites in the knownsite_pcdfs list. For now,
#'  only UTM coordinates are supported.
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
#'  It can be one of 'uniform' or 'exponential'—see 'weight' function.
#' @param weightparams A vector of parameters passed to the 'weight' function.
#' @param combinations A conveience parameter used to speed up calculations. The
#'  lamap calculation involves the use of the inclusion-exclusion principle for
#'  estimating the total probability of a union of independent events.
#'  Consequently, the R 'combn' function is used and can slow the calculations
#'  considerably. However, the user can instead use a utility function provided
#'  with this package—'prepCombinations'. It computes the combinations which can
#'  then be saved and passed to the lamap function with this parameter, saving
#'  the algorithm from having to regenerate the combinations for every cell in
#'  in the desired output map.
#' @return Returns '1' when finished executing. The LAMAP output surface will be
#'  stored in the 'outputpath'.
#' @export

parLamap <- function(cluster_object,
                     rasterpath,
                     outputpath,
                     knownsite_pcdfs,
                     knownsite_coords,
                     steps,
                     maxsites=NA,
                     weightfun,
                     weightparams,
                     combinations=NA,
                     nosupport=NA,
                     partial=T){
   rasterdata <- stack(rasterpath)
   lamap_surface <- raster(ext=extent(rasterdata),
                           crs=projection(rasterdata),
                           resolution=res(rasterdata))
   raster_output_cellnums <- matrix(1:ncell(lamap_surface),
                                    nrow=nrow(rasterdata),
                                    byrow=T)
   nrasterrows <- nrow(rasterdata)
   if(is.na(maxsites)){
      maxsites <- nrow(knownsite_coords)
   }
   l1 <- writeStart(lamap_surface,lamap_output_path,overwrite=T)
   for(j in 1:nrasterrows){
      lamaprow <- parSapply(cluster_object,
                        raster_output_cellnums[j,],
                        parLamapCaller,
                        rasterdata=rasterdata,
                        knownsite_pcdfs=knownsite_pcdfs,
                        knownsite_coords=knownsite_coords,
                        steps=steps,
                        maxsites=maxsites,
                        weightfun=weightfun,
                        weightparams=weightparams,
                        combinations=combinations,
                        nosupport=nosupport,
                        partial=partial)
      writeValues(l1,t(lamaprow),j)
   }
   l1 <- writeStop(l1)
   return(1)
}

#' parLamapCaller
#'
#' Wrapper designed to call the 'lamap' function for use in 'sapply'-like calls.
#'
#' @param x Int. Corresponds to the index of a single raster cell in a 'raster'
#' object. See the R package 'raster'.
#' @param rasterdata Raster object. Likely a raster stack from the package
#'  'raster'
#' @param knownsite_pcdfs List of of lists containing cdfs for each known site.
#'  See 'pcdf()'.
#' @param knownsite_coords Dataframe containing xy coordinates for each known
#'  site. The columns need to be named "id","x","y". The ids should be ordered
#'  to correspond to the order of sites in the knownsite_pcdfs list. For now,
#'  only UTM coordinates should be used.
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
#'  It can be one of 'uniform' or 'exponential'—see 'weight' function.
#' @param weightparams A vector of parameters passed to the 'weight' function.
#' @param combinations A conveience paramter used to speed up calculations. The
#'  lamap calculation involves the use of hte inclusion-exclusion principle for
#'  estimating the total probability of a union of independent events.
#'  Consequently, the R 'combn' function is used and can slow the calculations
#'  considerably. However, the user can instead use a utility function provided
#'  with this package—'prepCombinations'. It computes the combinations which can
#'  then be saved and passed to the lamap function with this parameter, saving
#'  the algorithm from having to regenerate the combinations for every cell in
#'  in the desired output map.
#' @return LAMAP value for one raster cell—see 'lamap'.
#' @export

parLamapCaller <- function(x,
                           rasterdata,
                           knownsite_pcdfs,
                           knownsite_coords,
                           steps,
                           maxsites,
                           weightfun,
                           weightparams,
                           combinations,
                           nosupport,
                           partial) {
   observedcell <- getObservedData(x, rasterdata)
   lamapoutput <- lamap(observed=observedcell,
                        knownsite_pcdfs=knownsite_pcdfs,
                        knownsite_coords=knownsite_coords,
                        steps=steps,
                        maxsites=maxsites,
                        weightfun=weightfun,
                        weightparams=weightparams,
                        combinations=combinations,
                        nosupport=nosupport,
                        partial=partial)
   return(lamapoutput)
}

#' getObservedData
#'
#' Extracts data for a single raster cell from a raster stack.
#'
#' @param x Cell index.
#' @param rasterdata A raster or stack object from the pacakge 'raster'.
#' @return A single row dataframe with columns corresponding to cell
#'  coordinates and variables.
#' @export

getObservedData <- function(x, rasterdata){
   observedcell <- as.data.frame(cbind(xyFromCell(rasterdata, x),
                                 round(extract(rasterdata, x))))
   return(observedcell)
}

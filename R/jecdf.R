#' multicdf
#'
#' Extract density estimates for observed data from multiple empirical
#' cumulative density (mass) functions. The function is designed to take one
#' set of observed values and return the corresponding density for each
#' variable.
#'
#' @param observed Single row dataframe of observed data (will be coerced).
#'  Each column should be one variable that corresponds to one of the cdfs in
#'  pcdfobj.
#' @param pcdfobj List of cdfs, each corresponding to a relevant predictive
#'  variable. See 'pcdf()'.
#' @param nosupport Value to return if the observed value(s) fall outside the
#'  support of the corresponding empirical cdf(s). Defaults to NA.
#' @param interpolate If True, the empirical cdf will be interpolated, providing
#'  an estimate of a continuous pdf instead. If False (default), then observed
#'  values that lie between mass density estimates will yield 'nosupport'.
#' @return Vector of densities (masses) corresponding to the variables in
#'  pcdfobj. The order of the densities (masses) will be the same as the order
#'  of the corresponding cdfs in pcdfobj.
#' @export

multicdf <- function(observed,pcdfobj,nosupport=NA,interpolate=F){
	mcdf_vals <- c()
	observed <- as.data.frame(observed)
	for (i in 1:ncol(observed)){
		training_densities_var <- pcdfobj$pcdfs[[i]]
		if(interpolate){
			min_training_value <- min(training_densities_var[,1])
			max_training_value <- max(training_densities_var[,1])
			training_densities_fun <- pcdfobj$cdfs_interpolated[[i]]
			insupport <- observed[,i] >= min_training_value & observed[,i] <= max_training_value
			if(insupport){
				cdf_value <- training_densities_fun(observed[,i])
			} else {
				cdf_value <- nosupport
			}
		} else {
			insupport <- observed[,i] %in% training_densities_var[,1]
			if(insupport){
				cdf_value <- training_densities_var[which(training_densities_var[,1]==observed[,i]),"cdfun"]
			} else {
				cdf_value <- nosupport
			}
		}
		mcdf_vals <- c(mcdf_vals,cdf_value)
	}
	return(mcdf_vals)
}

#' jdensity
#'
#' A wrapper for multicdf() that extracts densities (masses) for a given set
#' of observed values +/- the vector 'steps'.
#'
#' @param observed Single row dataframe of observed data (will be coerced).
#'  Each column should be one variable that corresponds to one of the cdfs in
#'  pcdfobj.
#' @param pcdfs List of cdfs, each corresponding to a relevant predictive
#'  variable. See 'pcdf()'.
#' @param nosupport Value to return if the observed value(s) fall outside the
#'  support of the corresponding empirical cdf(s). Defaults to NA.
#' @param partial If True, a density of >= 0 will be returned even if there is
#'  only partial overlap between the range (observed+/-steps) and the support
#'  of the corresponding cdf. The default value is True and it is generally
#'  recommended. It means that the function will return the density
#'  (mass) estimate for the portion of a given observed range that does fall
#'  within the support of the corresponding cdf.
#'  Without this option, the endpoint of the range that falls outside the
#'  cdf support will yield the nosupport value, which will mean the entire range
#'  could be disregraded in further calculations.
#' @param interpolate If True, the empirical cdf will be interpolated, providing
#'  an estimate of a continuous pdf instead. If False (default), then observed
#'  values that lie between mass density estimates will yield 'nosupport'.
#' @return Vector of densities (masses) corresponding to the variables in
#'  pcdfobj. The order of the densities (masses) will be the same as the order
#'  of the corresponding cdfs in pcdfobj.
#' @export

jdensity <- function(observed,steps,pcdfs,nosupport=NA,partial=T,interpolate=F){
	jdensity <- c()
	observed <- as.data.frame(observed)
	n_variables <- ncol(observed)
	observed_max <- observed + steps
	observed_min <- observed - steps
	pcdf_ranges <- lapply(pcdfs$pcdfs,function(x)range(x$value))
	pcdf_ranges <- do.call(cbind,pcdf_ranges)
	max_density <- multicdf(observed_max,pcdfs,nosupport,interpolate)
	min_density <- multicdf(observed_min,pcdfs,nosupport,interpolate)
	jdensities <- cbind(min_density,max_density)
	if(partial){
		for(j in 1:n_variables){
			overlap <- findOverlap(c(observed_min[j],observed_max[j]),pcdf_ranges[,j])
			if(overlap == 1){
				jdensities[j,1] <- 0
				jdensities[j,2] <- 1
			} else if(overlap == 2){
				jdensities[j,1] <- 0
			} else if(overlap == 3){
				jdensities[j,1] <- 0
				jdensities[j,2] <- 0
			} else if(overlap == 4){
				#leave jdensity alone
			} else if(overlap == 5){
				jdensities[j,2] <- 1
			} else if(overlap == 6){
				jdensities[j,1] <- 0
				jdensities[j,2] <- 0
			}
		}
	}
	return(jdensities)
}

#' pcdf
#'
#' Estimates the empirical probability and cumulative density (mass) functions
#' for a dataframe of raster cell values. This function is mostly used internally.
#'
#' @param traindf A dataframe of raster cells. Each row should be a single cell.
#'  Each column should be a variable of interest (e.g., elevation, slope, etc).
#' @param inteprolate Should the cdf(s) be interpolated to estimate a continuous
#'  pdf? Defaults to False. Ideally, all data should be rounded (effectively
#'  binned) before it's brought into R for lamap analysis. So, interpolation
#'  should not be used. It would undermine one of the fundamental advantages of
#'  of the lamap theory—namely, that multimodal terrain traits could be
#'  important for predictive modelling. Unless the data are very high resolution,
#'  interpolating could erroneously 'smooth out' the very features of the spatial
#'  cdfs that were important to site location.
#' @return List containing pcdfs.
#' @export

pcdf <- function(traindf,interpolate=F){
	traindf_dims <- dim(traindf)[1]
	traindf_densities <- densdf(traindf,traindf_dims[1])
	pdf_total_integrals <- unlist(lapply(traindf_densities,pdfIntegral))
	pcdfs <-  lapply(traindf_densities,function(x)data.frame(value=x[,1],pdfun=x[,2],cdfun=cumsum(x[,2])))
	if(interpolate){
		cdfs_interpolated <- lapply(pcdfs,function(x)splinefun(x$value,x$cdfun,method="hyman"))
	} else {
		cdfs_interpolated <- NA
	}
	return(list(pcdfs=pcdfs,total_integrals=pdf_total_integrals,cdfs_interpolated=cdfs_interpolated))
}

#' densdf
#'
#' Calculates density dataframes—used internally by other functions.
#'
#' @param df A dataframe of raster cells. Each row should be a single cell.
#'  Each column should be a variable of interest (e.g., elevation, slope, etc).
#' @return A list of dataframes containing the estimated empirical probability
#'  and cumulative density (mass) functions.
#' @export

densdf <- function(df,nobs){
   frequencies <- apply(df,2,function(x)as.data.frame(table(x)))
   frequencies <- lapply(frequencies,numericFactorLevels)
   densities <- lapply(frequencies,function(x)cbind(x[,1],x[,2]/nobs))
	densities <- lapply(densities,function(x){
		min_x <- min(x[,1])
		max_x <- max(x[,1])
		expanded_x <- c(min_x:max_x)
		densities_zero <- cbind(expanded_x,rep(0,length(expanded_x)))
		densities_zero[which(densities_zero[,1] %in% x[,1]),2]<-x[,2]
		data.frame(value=densities_zero[,1],dens=densities_zero[,2])
	})
   return(densities)
}

#' pdfIntegral
#'
#' Calculates a step-wise integral.
#'
#' @param df A dataframe with x and y columns. This function is used internally
#'  by the 'pcdf' function.
#' @return An integral estimate (single numeric value).
#' @export

pdfIntegral <- function(df){
   integral = sum(df[,2][1:(nrow(df)-1)] * diff(df[,1]))
   return(integral)
}

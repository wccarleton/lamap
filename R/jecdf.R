multicdf <- function(observed,pcdfobj,nosupport=NA,interpolate=F){
	mcdf_vals <- c()
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

jdensity <- function(observed,steps,pcdfs,nosupport=NA,partial=T,interpolate=F){
	jdensity <- c()
	observed_max <- observed + steps
	observed_min <- observed - steps
	max_density <- multicdf(observed_max,pcdfs,nosupport,interpolate)
	min_density <- multicdf(observed_min,pcdfs,nosupport,interpolate)
	jdensities <- cbind(min_density,max_density)
	if(partial){
		jdensities <- apply(jdensities,1,function(x){
			if(all(is.na(x))){
				return(x)
			} else if(is.na(x[1])){
				return(c(0,x[2]))
			} else if(is.na(x[2])){
				return(c(x[1],1))
			}
		})
	}
	return(jdensities)
}

pcdf <- function(traindf,interpolate=F){
	traindf_dims <- dim(traindf)[1]
	traindf_densities <- densdf(traindf,traindf_dims[1])
	print(traindf_densities)
	pdf_total_integrals <- unlist(lapply(traindf_densities,pdfIntegral))
	pcdfs <-  lapply(traindf_densities,function(x)data.frame(value=x[,1],pdfun=x[,2],cdfun=cumsum(x[,2])))
	if(interpolate){
		cdfs_interpolated <- lapply(pcdfs,function(x)splinefun(x$value,x$cdfun,method="hyman"))
	} else {
		cdfs_interpolated <- NA
	}
	return(list(pcdfs=pcdfs,total_integrals=pdf_total_integrals,cdfs_interpolated=cdfs_interpolated))
}

checkStepSizeEquality <- function(df){

}

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

pdfIntegral <- function(df){
   integral <- sum(df[,2][1:(nrow(df)-1)] * diff(df[,1]))
   return(integral)
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
         upper_index <- lower_index + 1
         lower_diff <- abs(x - v[lower_index])
         upper_diff <- abs(x - v[upper_index])
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

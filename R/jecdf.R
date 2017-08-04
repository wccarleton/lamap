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
			#print(overlap)
			#print(jdensities)
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
   integral = sum(df[,2][1:(nrow(df)-1)] * diff(df[,1]))
   return(integral)
}

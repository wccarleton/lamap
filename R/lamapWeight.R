weightKernel <- function(kernelType,ptsdist,param=c()) {
   maxdist = max(ptsdist['dist'])
	if (kernelType == "Uniform") {
		weights = ifelse(ptsdist['dist'] <= maxdist,1,0)
	}
	if (kernelType == "Triangle") {
		weights = ifelse(ptsdist['dist'] <= maxdist,1 - abs(ptsdist['dist']/maxdist),0)
	}
	if (kernelType == "Exponential") {
		if (length(param) < 1) {
			param[1] = 20 #dist declined to what weight
			#param[2] = 20 #rate of decline modifier
		}
		weights = exp((-ptsdist['dist']/maxdist)*param[1])
	}
	if (kernel == "Guassian") {
		if (length(param) > 0) {
			if (is.null(param[1])){param[1] = 1} #variance
			if (param[1] < 0.16){param[1] = 0.16} #variance limit to ensure upper weight limit is 1
			if (is.null(param[2])){param[1] = 1} #dispersion
			if (is.null(param[3])){param[1] = 3} #standard deviation
		}
		else {param=c(1,1,3)}
		x = (param[3] * sqrt(param[1])) / maxdist
		weights = (1 / sqrt(2 * pi * param[1])) * exp(-(((ptsdist['dist'] * x) / param[2])^2 / (2*param[1])))
	}
	colnames(weights) <- 'weight'
	return(cbind(ptsdist,weights))
}

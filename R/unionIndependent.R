#' Probability of union of independent events
#'
#' Calculates the probability of the union of independent events, given a vector of probabilities
#'
#' @param probs vector of probabilities
#'
#' @examples
#' unionIndependent(probs)
#'
#' @export
unionIndependent <- function(probs,combinations=NULL) {
	nprobs <- length(probs)
	intsctn <- c()
	if(!is.null(combinations)){
		for (i in 2:nprobs){
			outcomes <- apply(combinations[[i]],2,function(x,l)prod(l[x]),probs)
			intsctn[i] <- sum(outcomes)
		}
	} else{
		for (i in 2:nprobs){
			combinations <- combn(probs,i)
			outcomes <- apply(combinations,2,prod)
			intsctn[i] <- sum(outcomes)
		}
	}
	intsctn <- intsctn[2:(length(intsctn)-1)]
	n_intsctn <- length(intsctn)
	#the following calculation is based on the principle of inclusion-exclusion
	a <- sum(probs)
	b <- sum(intsctn[which(1:n_intsctn %% 2 == 1)])
	c <- sum(intsctn[which(1:n_intsctn %% 2 == 0)])
	d <- ((-1)^(n_intsctn+1) * prod(probs))
   prob_union <- a - b + c + d
	return(prob_union)
}

#' unionIndependent
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

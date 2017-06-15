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

unionIndependant <- function(probs) {
	intsctn <- c()
	column2list <- function(x){list(x)}
	for (i in 2:length(probs)){
		outcomes <- apply(combn(probs,i),2,column2list)
		for (j in 1:length(outcomes)){
				outcomes[[j]] <- outcomes[[j]][[1]]
		}
		outcomes_container <- unlist(lapply(outcomes,prod))#unlist(mpi.parLapply(outcomes,prod))
		intsctn[i] <- sum(outcomes_container)
	}
	intsctn <- intsctn[-1]
   prob_union <- sum(probs) - sum(intsctn[which(1:length(intsctn) %% 2 == 1)]) + sum(intsctn[which(1:length(intsctn) %% 2 == 0)]) + ((-1)^length(intsctn) * prod(probs))
	return(prob_union)
}

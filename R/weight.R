#' weight
#'
#' Weighting function for the 'lamap' algorithm. This function is essentially a
#' switch that passes 'params' and 'distances' to the appropriate weighting
#' function, which then returns a vector with the same length as 'distances'
#' containing the weights.
#'
#' Uniform: Uniform weighting based on distance. The 'param' vector should be a
#' a single value specifying the maximum distance (in meters), beyond which
#' weight is 0.
#'
#' Exponential: Exponential weighting where the first element in 'param' affects
#' the rate of decline in weight with increasing distance—higher values increase
#' rate. The second element specifies the distance at which the weight is very
#' nearly zero. The default for the second element is the maximum distance in
#' the vector 'distances'. The exponential function is as follows:
#' exp((-x/param[2])*param[1])
#'
#' @param distances Vector of distances—returned from 'orderSites'.
#' @param weightfun The name of the weight function to use. Currently can be one
#'  of 'uniform' or 'exponential'
#' @param param A vector of parameters to be used in the weighting function. See
#'  the apropriate function for details.
#' @return Vector of weights in the same order as 'distances'.
#' @export

weight <- function(distances,weightfun,param=c()){
   switch(weightfun,
      uniform = w_uniform(distances,param),
      exponential = w_exponential(distances,param))
}

#' w_uniform
#'
#' Calculates a vector of weights using a uniform weighting function.
#'
#' @param x A vector of distances.
#' @param param A vector of weighting function parameters. See primary 'weight'
#' function.
#' @return A length n list of matrices.
#' @export

w_uniform <- function(x,param){
   y <- ifelse(x <= param[1],1,0)
   return(y)
}

#' w_exponential
#'
#' Calculates a vector of weights using an exponential weighting function.
#'
#' @param x A vector of distances.
#' @param param A vector of weighting function parameters. See primary 'weight'
#' function.
#' @return A length n list of matrices.
#' @export

w_exponential <- function(x,param){
   if(!is.na(param[2]) || !is.null(param[2])){
      max_x <- param[2]
   } else {
      max_x <- max(x)
   }
   y = exp((-x/max_x)*param[1])
   return(y)
}

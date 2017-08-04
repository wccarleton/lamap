#' LAMAP weighting functions
#'
#'
#'
#'
#'
#'
weight <- function(distances,weightfun,param=c()){
   switch(weightfun,
      uniform = w_uniform(distances,param),
      exponential = w_exponential(distances,param))
}

w_uniform <- function(x,param){
   y <- ifelse(x <= param[1],1,0)
   return(y)
}

w_exponential <- function(x,param){
   if(!is.na(param[2])){
      max_x <- param[2]
   } else {
      max_x <- max(x)
   }
   y = exp((-x/max_x)*param[1])
   return(y)
}

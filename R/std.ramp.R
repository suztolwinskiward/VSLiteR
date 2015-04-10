#' "Standard ramp" function for building growth response functions.
#' 
#' \code{std.ramp.r} takes a lower and upper bounds, and creates a linear response
#' between 0 and 1 for values in between the bounds.  Values below (above) the lower (upper) 
#' bound are assigned a value of zero (one).
#' 
#' @param x The value at which we want to evaluate the ramp function.
#' @param x1 The lower bound of the support of the nonzero part of the ramp.
#' @param x2 The lower bound of the range of the preimage of 1.
#' 
#' @export

std.ramp <- function(x,x1,x2){return(
  apply(as.matrix(apply((x-x1)/(x2-x1),1:length(dim(x)),min,1)),
        1:length(dim(x)),max,0)
)}



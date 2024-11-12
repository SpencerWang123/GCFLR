#' Huber loss function to transform residuals
#'
#' @param x residual term of term model
#' @param c huber parameter
#'
#' @return transformed values

huber_loss=function(x,c=1.345)
  {return(0.5*x^2*(abs(x)<=c)+(c*abs(x)-0.5*c^2)*(abs(x)>c))}

#' multinomial parameters for the full condition distribution at the mask position
#'
#' @param reg_dt N*M principle components scores of the observations
#' @param estimated_parameters coefficient estimation of the clustering regression results
#' @param mask_index int in 1:N used to mask several specific label (for full condition distribution estimation)
#' @param K number of clusters
#' @param lambda parameters used to construct a multinomial distribution
#' @param method loss function used in the regression
#'
#' @return multinomial parameters for the full condition distribution at the mask position

distribution_parameters=function(reg_dt,estimated_parameters,mask_index,K,lambda=1,method="ls")
{
  temp_observation=c(1,as.numeric(reg_dt[mask_index,-1]))
  temp_y=reg_dt[mask_index,1]
  y_estimated=apply(temp_observation*estimated_parameters,2,sum)
  fun_745=function(x){return(min(x,200))}
  if(method=="ls")
  {parameters_multinomial=exp(-1*lambda*sapply((y_estimated-temp_y)^2,fun_745))/sum(exp(-1*lambda*sapply((y_estimated-temp_y)^2,fun_745)))}
  if(method=="hl")
  {parameters_multinomial=exp(-1*lambda*huber_loss(y_estimated-temp_y))/sum(exp(-1*lambda*huber_loss(y_estimated-temp_y)))}
  if(method=="lad")
  {parameters_multinomial=exp(-1*lambda*abs(y_estimated-temp_y))/sum(exp(-1*lambda*abs(y_estimated-temp_y)))}
  if(method=="logit")
  {parameters_multinomial=exp(-1*lambda*sapply((y_estimated-temp_y)^2,fun_745))/sum(exp(-1*lambda*sapply((y_estimated-temp_y)^2,fun_745)))}

  return(parameters_multinomial)
}


#' BIC loss result to determine the number of clusters
#'
#' @param reg_dt N*M principle components scores of the observations
#' @param label N*1 label vector that conveys clustering results
#' @param K number of clusters
#' @param method loss function used in the regression
#'
#' @return calculated BIC results

BIC_fun=function(reg_dt,label,temp_K_seq,method="ls")
{
  p=ncol(reg_dt)#number of regression variables including constant variable(intercept)
  dt_list=NULL;for(i in temp_K_seq){dt_list=c(dt_list,list(reg_dt[label==i,]))}
  #BIC penalty term
  n=nrow(reg_dt)
  num_par=ncol(reg_dt)
  qk=max(temp_K_seq)*ncol(reg_dt)#K*P
  An=log(log(n))/n
  penalty_term=qk*An
  # LS
  if(method=="ls")
  {residuals_reg=unlist(lapply(dt_list, ls_reg_loss,p=p))
  BIC_res=mean(residuals_reg^2)+5*penalty_term
  }
  # HL
  if(method=="hl")
  {residuals_reg=unlist(lapply(dt_list, hl_reg_loss,p=p))
  BIC_res=log(mean(huber_loss(residuals_reg)))+5*penalty_term
  }
  # LAD
  if(method=="lad")
  {
    residuals_reg=unlist(lapply(dt_list, lad_reg_loss,p=p))
    BIC_res=log(mean(abs(residuals_reg)))+5*penalty_term
  }
  # Logit
  if(method=="logit")
  {
    residuals_reg=unlist(lapply(dt_list, logit_reg_loss,p=p))
    BIC_res=log(mean(abs(residuals_reg)))+5*penalty_term
  }
  return(BIC_res)#first row of the matrix is intercept
}

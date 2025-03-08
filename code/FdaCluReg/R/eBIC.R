#' eBIC loss result to determine the number of clusters
#'
#' @param reg_dt N*M principle components scores of the observations
#' @param label N*1 label vector that conveys clustering results
#' @param K number of clusters
#' @param method loss function used in the regression
#'
#' @return calculated BIC results

eBIC_fun=function(reg_dt,label,temp_K_seq,method="ls",xi_1=1,xi_2=0.5,xi_3=1)
{
  library(gmp)
  p=ncol(reg_dt)#number of regression variables including constant variable(intercept)
  dt_list=NULL;for(i in temp_K_seq){dt_list=c(dt_list,list(reg_dt[label==i,]))}#dt_list是list按照类簇编码进行排列
  #BIC penalty term
  n=nrow(reg_dt)
  num_par=ncol(reg_dt)
  k_num=max(temp_K_seq)
  qk=k_num*num_par#K*P
  An=log(n)
  logS=log(Stirling2(n,k_num))

  penalty_term=(xi_1*logS+xi_2*qk*(An^xi_3))
  # LS
  if(method=="ls")
  {residuals_reg=unlist(lapply(dt_list, ls_reg_loss,p=p))
  eBIC_res=mean(residuals_reg^2)+penalty_term/n
  }
  # HL
  if(method=="hl")
  {residuals_reg=unlist(lapply(dt_list, hl_reg_loss,p=p))
  eBIC_res=log(mean(huber_loss(residuals_reg)))+penalty_term/n
  }
  # LAD
  if(method=="lad")
  {
    residuals_reg=unlist(lapply(dt_list, lad_reg_loss,p=p))
    eBIC_res=log(mean(abs(residuals_reg)))+penalty_term/n
  }
  # Logit
  if(method=="logit")
  {
    residuals_reg=unlist(lapply(dt_list, logit_reg_loss,p=p))
    eBIC_res=log(mean(abs(residuals_reg)))+penalty_term/n
  }

  return(eBIC_res)#first row of the matrix is intercept
}

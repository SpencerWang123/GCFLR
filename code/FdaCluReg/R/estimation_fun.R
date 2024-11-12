#' coefficient estimation of the clustering regression results
#'
#' @param reg_dt N*M principle components scores of the observations
#' @param label N*1 label vector that conveys clustering results
#' @param mask_index int in 1:N used to mask several specific label (for full condition distribution estimation)
#' @param method loss function used in the regression
#' @param K number of clusters
#'
#' @return list(parameters_mat, residuals_list):
#' parameters_mat: (M+1)*K coefficient matrix whose first row indicates the intercept term
#' residuals_list: residuals of the estimation

estimation_fun=function(reg_dt,label,mask_index,method="ls",temp_K_seq)
{
  p=ncol(reg_dt)#number of regression variables including constant variable(intercept)
  if(!is.null(mask_index))
  {reg_dt_temp=reg_dt[-mask_index,];label_temp=label[-mask_index]}else{reg_dt_temp=reg_dt;label_temp=label}

  dt_list=NULL;for(i in temp_K_seq){dt_list=c(dt_list,list(reg_dt_temp[label_temp==i,]))}
  #coefficient estimation
  if(method=="ls"){parameters=lapply(dt_list, ls_reg_fun,p=p)}
  if(method=="hl"){parameters=lapply(dt_list, hl_reg_fun,p=p)}
  if(method=="lad"){parameters=lapply(dt_list, lad_reg_fun,p=p)}
  if(method=="logit"){parameters=lapply(dt_list, logit_reg_fun,p=p)}
  parameters_mat=matrix(unlist(parameters),ncol = length(temp_K_seq))
  colnames(parameters_mat)=temp_K_seq
  #residual outputs
  if(method=="ls"){residuals_res=lapply(dt_list, ls_reg_loss,p=p)}
  if(method=="hl"){residuals_res=lapply(dt_list, hl_reg_loss,p=p)}
  if(method=="lad"){residuals_res=lapply(dt_list, lad_reg_loss,p=p)}
  if(method=="logit"){residuals_res=lapply(dt_list, logit_reg_loss,p=p)}
  return(list(parameters_mat=parameters_mat,residuals_list=residuals_res))#first row of the matrix is intercept
}

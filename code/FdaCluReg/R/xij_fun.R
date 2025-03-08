#' generating simulated X_{ij}(i=1,2,...,N)
#'
#' @param function_list K function objects used in simulation compiled into a list
#' @param fda_coef_df coefficient matrix with K columns and N rows
#' @param sample_t T*1 vector indicating observing time points
#' @param error_sd scalar indicating variance of normal observing error
#' @param seed random seed
#'
#' @return simulated N*T matrix of sample X_i at time t_j
#'
xij_fun=function(function_list,fda_coef_df,sample_t,error_sd,seed=1)
{
  set.seed(seed)
  #基函数值
  phi_df=NULL
  for (i in 1:length(function_list) )#计算phi(tj)
  {temp_fun=function_list[i][[1]];phi_df=c(phi_df,temp_fun(sample_t))}
  phi_df=t(matrix(phi_df,ncol = length(function_list)))#K个变量，K*T的矩阵

  #组合样本,包含偏差
  xij_mat=NULL
  for(j in 1:nrow(fda_coef_df))
  {xij=apply(fda_coef_df[j,]*phi_df,2,sum);xij=xij+rnorm(length(xij),0,error_sd);xij_mat=c(xij_mat,xij)}
  xij_mat=matrix(xij_mat,ncol = nrow(fda_coef_df))
  return(xij_mat)
}

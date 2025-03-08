#' generating R rounds of sample results of labels of the N observations
#'
#' @param reg_dt N*M principle components scores of the observations
#' @param label_initial N*1 vector of initialized label sequence
#' @param lambda scalar indicating parameters used to construct a multinomial distribution
#' @param round_num R rounds-scalar indicating how many round the algorithm runs
#' @param seed random seed
#' @param method loss function used in the regression
#'
#' @return N*R label matrix indicating the sampling results of each round

gibbs_sample=function(reg_dt,label_initial,lambda=0.01,round_num=1000,seed=1,method="ls")
{
  set.seed(seed)
  label_list=NULL#vector to  save samples
  K=max(label_initial)# number of groups
  N=nrow(reg_dt)# number of samples
  temp_label=label_initial
  for(r in 1:round_num)
  {
    for (temp_mask_index in 1:N)
    {
      temp_parameters=estimation_fun(reg_dt = reg_dt,label = temp_label,mask_index = temp_mask_index,method=method,temp_K_seq = 1:K)$parameters_mat
      temp_multinomial=distribution_parameters(reg_dt = reg_dt,estimated_parameters = temp_parameters,mask_index = temp_mask_index,lambda=lambda,method = method,K=K)
      sample_label=as.numeric(rmultinom(1,1,temp_multinomial))
      temp_label[temp_mask_index]=which(sample_label==1)
      #print(temp_parameters);print(temp_multinomial)
    }
    label_list=c(label_list,temp_label)
  }

  label_mat=matrix(label_list,ncol = round_num)
  return(label_mat)
}

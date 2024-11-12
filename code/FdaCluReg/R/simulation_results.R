simulation_results_function=function(X_list,Y_list,basis,sample_t,var_ratio=0.9,
                                     K_max=10,true_label=NULL,rounds=110,burning_period=10,lambda=1,method,ncomp_u=10,ncomp_l=1,seed,type=c(TRUE,TRUE))
{
  res_gibbs=NULL
  res_2021=NULL
  assement_all_BIC=NULL
  assement_all_eBIC=NULL
  selected_K_gibbs_BIC=NULL
  selected_K_2021_BIC=NULL
  selected_K_gibbs_eBIC=NULL
  selected_K_2021_eBIC=NULL
  N=length(Y_list[[1]])
  if(is.null(true_label)){true_label=rep(1,N)}
  for (i in 1:length(Y_list))
  {
    if(type[1])
    {
      temp_gibbs=gibbs_clusterwise_regression_auto(X_list=X_list,
                                                   Y=Y_list[[i]],
                                                   basis=basis,
                                                   sample_t=sample_t,
                                                   var_ratio=var_ratio,
                                                   K_max=K_max,
                                                   true_label=true_label,
                                                   rounds=rounds,
                                                   burning_period=burning_period,
                                                   lambda = lambda,
                                                   method=method,
                                                   ncomp_u=10,ncomp_l=ncomp_l,seed)
      res_gibbs=c(res_gibbs,list(temp_gibbs))
      selected_K_gibbs_BIC=c(selected_K_gibbs_BIC,temp_gibbs$K_BIC)
      selected_K_gibbs_eBIC=c(selected_K_gibbs_eBIC,temp_gibbs$K_eBIC)
    }  #saving

    if(type[2])
    {
      temp_2021=clusterwise_reg_2021_auto(X_list =X_list,
                                          Y=Y_list[[i]],
                                          basis=basis,
                                          sample_t=sample_t,
                                          var_ratio=var_ratio,
                                          K_max=K_max,
                                          true_label=true_label,
                                          ncomp_u=10,ncomp_l=ncomp_l,seed,method = method)
      res_2021=c(res_2021,list(temp_2021))
      selected_K_2021_BIC=c(selected_K_2021_BIC,temp_2021$K_BIC)
      selected_K_2021_eBIC=c(selected_K_2021_eBIC,temp_2021$K_eBIC)
      }  #saving

    if(!type[1]){temp_gibbs=temp_2021}
    assement_all_BIC=c(assement_all_BIC,temp_gibbs$selected_res_BIC$assessment,temp_2021$selected_res_BIC$assessment)
    assement_all_eBIC=c(assement_all_eBIC,temp_gibbs$selected_res_eBIC$assessment,temp_2021$selected_res_eBIC$assessment)
  }
  if(type[1]){names(res_gibbs)=c("homo_n","homo_t","homo_chisq","heter_n","heter_t","heter_chisq")}
  if(type[2]){names(res_2021)=c("homo_n","homo_t","homo_chisq","heter_n","heter_t","heter_chisq")}
  #assement
  assement_all_mat_BIC=matrix(assement_all_BIC,byrow = T,ncol = 3);assement_all_mat_eBIC=matrix(assement_all_eBIC,byrow = T,ncol = 3)
  colnames(assement_all_mat_BIC)=c("ARI","RI","NMI");colnames(assement_all_mat_eBIC)=c("ARI","RI","NMI")
  rownames(assement_all_mat_BIC)=1:nrow(assement_all_mat_BIC);rownames(assement_all_mat_eBIC)=1:nrow(assement_all_mat_eBIC)

  rownames(assement_all_mat_BIC)[(1:6)*2-1]=paste("gibbs",c("homo_n","homo_t","homo_chisq","heter_n","heter_t","heter_chisq"),sep = "_")
  rownames(assement_all_mat_BIC)[(1:6)*2]=paste("m2021",c("homo_n","homo_t","homo_chisq","heter_n","heter_t","heter_chisq"),sep = "_")

  rownames(assement_all_mat_eBIC)[(1:6)*2-1]=paste("gibbs",c("homo_n","homo_t","homo_chisq","heter_n","heter_t","heter_chisq"),sep = "_")
  rownames(assement_all_mat_eBIC)[(1:6)*2]=paste("m2021",c("homo_n","homo_t","homo_chisq","heter_n","heter_t","heter_chisq"),sep = "_")

  res_all=list(res_gibbs_all=res_gibbs,res_2021_all=res_2021,
               assement_all_mat_BIC=assement_all_mat_BIC,
               assement_all_mat_eBIC=assement_all_mat_eBIC,
               selected_K_2021_BIC=c(selected_K_gibbs_BIC,selected_K_2021_BIC),
               selected_K_2021_eBIC=c(selected_K_gibbs_eBIC,selected_K_2021_eBIC))
  return(res_all)
}

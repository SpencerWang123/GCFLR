gibbs_clusterwise_regression_auto=function(X_list,Y,basis, sample_t,var_ratio = 0.9,K_max,true_label,
                                           rounds = 50, burning_period = 10,method = "ls",lambda = 0.1,
                                           ncomp_u = 10,ncomp_l = 2,seed = 1)
{
  K_seq=2:K_max
  test_res=gibbs_clusterwise_regression(X_list=X_list,Y=Y,basis=basis,sample_t=sample_t,var_ratio=var_ratio,
                                        K_seq=K_seq,true_label=true_label,rounds=rounds,
                                        burning_period=burning_period,lambda=lambda,method=method,
                                        ncomp_u=ncomp_u,ncomp_l=ncomp_l,seed=seed)
  selected_K_BIC=K_seq[which.min(test_res$BIC)]
  selected_K_eBIC=K_seq[which.min(test_res$eBIC)]
  selected_res_BIC=gibbs_clusterwise_regression(X_list=X_list,Y=Y,basis=basis,sample_t=sample_t,var_ratio=var_ratio,
                                            K_seq=selected_K_BIC,true_label=true_label,rounds=rounds,
                                            burning_period=burning_period,lambda=lambda,method=method,
                                            ncomp_u=ncomp_u,ncomp_l=ncomp_l,seed=seed)
  selected_res_eBIC=gibbs_clusterwise_regression(X_list=X_list,Y=Y,basis=basis,sample_t=sample_t,var_ratio=var_ratio,
                                            K_seq=selected_K_eBIC,true_label=true_label,rounds=rounds,
                                            burning_period=burning_period,lambda=lambda,method=method,
                                            ncomp_u=ncomp_u,ncomp_l=ncomp_l,seed=seed)
  return(list(selected_res_BIC=selected_res_BIC,K_BIC=selected_K_BIC,
              selected_res_eBIC=selected_res_eBIC,K_eBIC=selected_K_eBIC,all_res=test_res))
}


clusterwise_reg_2021_auto=function(X_list,Y,basis,sample_t,var_ratio,K_max,true_label,
                                   method = "ls",ncomp_u = 10,ncomp_l = 2,seed = 1)
{
  K_seq=2:K_max
  test_res=clusterwise_reg_2021(X_list=X_list,Y=Y,basis=basis,sample_t=sample_t,var_ratio=var_ratio,K_seq=K_seq,
                                true_label=true_label,method=method,ncomp_u=ncomp_u,ncomp_l=ncomp_l,seed=seed)
  selected_K_BIC=K_seq[which.min(test_res$BIC)]
  selected_K_eBIC=K_seq[which.min(test_res$eBIC)]
  selected_res_BIC=clusterwise_reg_2021(X_list=X_list,Y=Y,basis=basis,sample_t=sample_t,var_ratio=var_ratio,K_seq=selected_K_BIC,
                                    true_label=true_label,method=method,ncomp_u=ncomp_u,ncomp_l=ncomp_l,seed=seed)
  selected_res_eBIC=clusterwise_reg_2021(X_list=X_list,Y=Y,basis=basis,sample_t=sample_t,var_ratio=var_ratio,K_seq=selected_K_eBIC,
                                    true_label=true_label,method=method,ncomp_u=ncomp_u,ncomp_l=ncomp_l,seed=seed)
  return(list(selected_res_BIC=selected_res_BIC,K_BIC=selected_K_BIC,
              selected_res_eBIC=selected_res_eBIC,K_eBIC=selected_K_eBIC,all_res=test_res))
}

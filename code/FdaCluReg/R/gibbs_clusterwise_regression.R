#' integrated gibbs sampling and estimation algorithm for clustering regression
#'
#' @param X_list list of length K. each element is a N*T matrix. contains simulated X_1:K observed at T evenly spaced time points
#' @param Y N*1 vector of the N response variable of the N observations
#' @param basis basis used in the smoothing
#' @param sample_t T evenly spaced time points of the observations
#' @param var_ratio minimum variance demanded to determine how many components are preserved in the FPCA
#' @param K_seq int. vector contains the numbers of clusters
#' @param true_label N*1 label vector contains ground truth labels
#' @param rounds R rounds-scalar indicating how many round the algorithm runs
#' @param burning_period scalar indicating how many rounds of results are abandoned in the sampling process
#' @param method loss function used in the regression with candidates c("ls","hl","lad")
#' @param lambda scalar indicating the parameters used to construct a multinomial distribution
#' @param ncomp_u maximum number of components demanded in FPCA
#' @param ncomp_l minimum number of components demanded in FPCA
#' @param seed random seed
#'
#' @return list(assessment of the results,BIC_seq,coefficients_list,labels_trace,fpca,components scores,beta_trace,residuals_trace,residuals_list,beta_fd)
#' assessment of the results: ARI and other assessment standards.
#' BIC_seq: sequence of BIC with the same length as K_seq
#' coefficients_list: list consisting of the coefficients matrix estimated under each number of clusters with the same length as K_seq.
#' label trace: trace of the label sequence throughout sampling process
#' fpca: fuinctional principle analysis results of X_i(t_j)
#' components scores: principle components scores of each observations
#' beta trace: each label sequence of the label trace generates one coefficient estimation
#' residual trace: each label sequence of the label trace generates one residuals sequence.
#' beta fd: estimated coefficient in the form of fd object

gibbs_clusterwise_regression=function(X_list,
                                      Y,
                                      basis,sample_t,var_ratio=0.9,
                                      K_seq,
                                      true_label,
                                      rounds=50,
                                      burning_period=10,
                                      method="ls",
                                      lambda=0.1,
                                      ncomp_u=10,ncomp_l=2,seed=1,need_trace=F)
{## integrated algorithm
  set.seed(seed)
  library(fda)
  # FDA
  ## smoothing variable by variable
  N=length(Y)
  pc_scores_all=data.frame(index=1:N)#empty dataframe for saving
  pca_res=NULL
  for (variable_index in 1:length(X_list))
  {
    X=X_list[variable_index][[1]]
    smoothed_x=smooth.basis(sample_t,X,basis)#fd object of X_ij
    ## fpca
    pca_x=pca.fd(smoothed_x$fd,nharm = ncomp_u,centerfns = F)# implement of FPCA
    M=which(cumsum(pca_x$varprop)>0.9)[1];if(M<ncomp_l){M=ncomp_l}# choose number of principle components
    pc_scores=pca_x$scores[,1:M]
    pca_res=c(pca_res,list(pca_x))# save fpca
    pc_scores_all=data.frame(pc_scores_all,pc_scores)
  }
  pc_scores_all=pc_scores_all[,-1]#delete index column

  # sampling process
  # results saving
  assess_df=NULL
  BIC_seq=NULL
  eBIC_seq=NULL
  label_df=NULL
  label_trace=NULL
  coefficient_list=NULL
  beta_trace=NULL
  residuals_trace=NULL
  residuals_list=NULL
  beta_fd_list=NULL
  ## for different number of groups
  for (K in K_seq)
  {
    init_label=sample(1:K,N,replace = T)# random initialization
    reg_dt=data.frame(y=Y,pc_scores_all)# data for regression
    #sample process
    if(method=="ls"){sample_res=gibbs_sample(reg_dt = reg_dt,label_initial = init_label,round_num = rounds,method = "ls",lambda = lambda)}
    if(method=="lad"){sample_res=gibbs_sample(reg_dt = reg_dt,label_initial = init_label,round_num = rounds,method = "lad",lambda = lambda)}
    if(method=="hl"){sample_res=gibbs_sample(reg_dt = reg_dt,label_initial = init_label,round_num = rounds,method = "hl",lambda = lambda)}
    if(method=="logit"){sample_res=gibbs_sample(reg_dt = reg_dt,label_initial = init_label,round_num = rounds,method = "logit",lambda = lambda)}
    sample_res_plot=sample_res
    #label_trace
    label_trace=c(label_trace,list(sample_res_plot))
    #estimation_trace
    if(need_trace)
    {
    estimation_trace=apply(sample_res_plot,2,estimation_fun,reg_dt=reg_dt,temp_K_seq=1:K,mask_index=NULL)
    #beta trace
    beta_trace_mat=vector("list",rounds)
    for (i in 1:ncol(sample_res_plot))
    {beta_trace_mat[[i]]=estimation_trace[[i]][[1]]}
    names(beta_trace_mat)=paste("round",1:rounds)
    beta_trace=c(beta_trace,list(beta_trace_mat))

    #residuals trace
    residuals_trace_mat=vector("list",rounds)
    for (i in 1:ncol(sample_res_plot))
    {residuals_trace_mat[[i]]=estimation_trace[[i]][[2]]}
    names(residuals_trace_mat)=paste("round",1:rounds)
    residuals_trace=c(residuals_trace,list(residuals_trace_mat))
    }else{
      beta_trace=c(beta_trace,list("no record"))
      residuals_trace=c(residuals_trace,list("no record"))
    }


    #burning period
    preserved_rounds=burning_period:ncol(sample_res)
    mode_group=apply(sample_res[,preserved_rounds],1,find_mode)

    # estimation
    estimation_temp=estimation_fun(reg_dt,mode_group,mask_index=NULL,method=method,temp_K_seq=1:K)
    coef_temp=estimation_temp$parameters_mat
    residuals_temp=estimation_temp$residuals_list

    #fd coefficient
    coef_fd=matrix(as.vector(coef_temp[-1,]),nrow = M)
    coef_fd_transform=t(t(coef_fd)%*%t(pca_x$harmonics$coefs[,1:M]))
    temp_beta_fd=fd(coef_fd_transform,pca_x$harmonics$basis)
    varibale_num=ncol(coef_fd_transform)/K
    temp_fd_list=vector("list",K)
    for (clu in 1:K){temp_fd_list[[clu]]=temp_beta_fd[(clu-1)*varibale_num+1:varibale_num]}
    names(temp_fd_list)=paste("cluster",1:K,sep = "-")
    beta_fd_list=c(beta_fd_list,list(temp_fd_list))

    #最终标签结果判定-贪婪搜索
    estimated_y=as.matrix(cbind(rep(1,N),reg_dt[,2:ncol(reg_dt)]))%*%coef_temp
    abs_error=abs(estimated_y-reg_dt[,1])
    final_label=apply(abs_error,1,which.min)


    # results assessment
    assess_res=c(ari(final_label,true_label),ri(final_label,true_label),nmi(final_label,true_label))

    #BIC
    BIC_res=BIC_fun(reg_dt = reg_dt,final_label,temp_K_seq = 1:K,method = method)
    eBIC_res=eBIC_fun(reg_dt = reg_dt,final_label,temp_K_seq = 1:K,method = method)

    #saving results
    assess_df=c(assess_df,assess_res)
    BIC_seq=c(BIC_seq,BIC_res)
    eBIC_seq=c(eBIC_seq,eBIC_res)
    label_df=c(label_df,final_label)
    coefficient_list=c(coefficient_list,list(coef_temp))
    residuals_list=c(residuals_list,list(residuals_temp))
  }

  # reshape results
  ## assess_df
  assess_df=data.frame(t(matrix(assess_df,ncol = length(K_seq))));
  rownames(assess_df)=paste(method,K_seq,sep = "-")
  colnames(assess_df)=c("ARI","RI","NMI")
  ## BIC_seq
  names(BIC_seq)=paste("cluster number=",K_seq)
  names(eBIC_seq)=paste("cluster number=",K_seq)
  ## label_df
  label_df=data.frame(matrix(label_df,ncol = length(K_seq)))
  rownames(label_df)=1:N
  colnames(label_df)=paste("group_result",K_seq,sep = "-")
  ## coefficient_list
  names(coefficient_list)=paste("cluster number=",K_seq)
  ## residuals_list
  names(residuals_list)=paste("cluster number=",K_seq)
  ##beta trace
  names(beta_trace)=paste("cluster number=",K_seq)
  ##residuals trace
  names(residuals_trace)=paste("cluster number=",K_seq)
  ##beta_fd_list
  names(beta_fd_list)=paste("cluster number=",K_seq)
  ##label_trace
  names(label_trace)=paste("cluster number=",K_seq)
  # outputs
  results_algorithm=list(assessment=assess_df,BIC=BIC_seq,eBIC=eBIC_seq,coefficients=coefficient_list,
                         labels=label_df,label_trace=label_trace,fpca=pca_res,reg_dt=reg_dt,
                         beta_trace=beta_trace,residuals_trace=residuals_trace,
                         residuals_list=residuals_list,beta_fd=beta_fd_list)
  cluster_selection=cluster_number_selection(results_algorithm)
  return(list(assessment=assess_df,BIC=BIC_seq,eBIC=eBIC_seq,coefficients=coefficient_list,labels=label_df,label_trace=label_trace,fpca=pca_res,reg_dt=reg_dt,
              beta_trace=beta_trace,residuals_trace=residuals_trace,residuals_list=residuals_list,beta_fd=beta_fd_list,cluster_selection=cluster_selection))
}

#' clusterwise regression algorithm proposed in 2021
#'
#' @param X_list list of length K. each element is a N*T matrix. contains simulated X_1:K observed at T evenly spaced time points
#' @param Y N*1 vector of the N response variable of the N observations
#' @param basis basis used in the smoothing
#' @param sample_t T evenly spaced time points of the observations
#' @param var_ratio minimum variance demanded to determine how many components are preserved in the FPCA
#' @param K_seq int. vector contains the numbers of clusters
#' @param true_label N*1 label vector contains ground truth labels
#' @param method loss function used in the regression with candidates c("ls","hl","lad")
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


clusterwise_reg_2021=function(X_list,
                              Y,
                              basis,sample_t,var_ratio,
                              K_seq,
                              true_label,
                              method="ls",
                              ncomp_u=10,ncomp_l=2,seed=1)
{
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
  N=nrow(as.matrix(pc_scores))

  # results saving
  assess_df=NULL
  BIC_seq=NULL
  eBIC_seq=NULL
  label_df=NULL
  coefficient_list=NULL
  beta_trace=NULL
  residuals_trace=NULL
  residuals_list=NULL
  beta_fd_list=NULL
  label_trace=NULL

  ## for different number of groups
  for (K in K_seq)
  {
    init_label=sample(1:K,N,replace = T)# random initialization
    reg_dt=data.frame(y=Y,pc_scores_all)# data for regression
    results=clustering_fun2021(reg_dt=reg_dt,label=init_label,method=method,K=K)
    #label_trace
    label_trace=c(label_trace,list(results$label_record))


    #beta trace
    beta_trace=c(beta_trace,list(results$beta_trace))

    #residuals trace
    residuals_trace=c(residuals_trace,list(results$residuals_trace))

    # coefficient estimation
    coef_temp=results$cofficients
    residuals_temp=results$residuals_temp
    names(residuals_temp)=paste("cluster",1:K,sep = "-")
    residuals_list=c(residuals_list,list(residuals_temp))


    #fd coefficient
    coef_fd=matrix(as.vector(coef_temp[-1,]),nrow = M)
    coef_fd_transform=t(t(coef_fd)%*%t(pca_x$harmonics$coefs[,1:M]))
    temp_beta_fd=fd(coef_fd_transform,pca_x$harmonics$basis)
    varibale_num=ncol(coef_fd_transform)/K
    temp_fd_list=vector("list",K)
    for (clu in 1:K){temp_fd_list[[clu]]=temp_beta_fd[(clu-1)*varibale_num+1:varibale_num]}
    names(temp_fd_list)=paste("cluster",1:K,sep = "-")
    beta_fd_list=c(beta_fd_list,list(temp_fd_list))


    # results assessment
    label_est=results$label
    assess_res=c(ari(label_est,true_label),ri(label_est,true_label),nmi(label_est,true_label))

    #BIC
    BIC_res=BIC_fun(reg_dt = reg_dt,label_est,temp_K_seq = 1:K,method = method)
    eBIC_res=eBIC_fun(reg_dt = reg_dt,label_est,temp_K_seq = 1:K,method = method)


    #saving results
    assess_df=c(assess_df,assess_res)
    BIC_seq=c(BIC_seq,BIC_res)
    eBIC_seq=c(eBIC_seq,eBIC_res)
    label_df=c(label_df,label_est)
    coefficient_list=c(coefficient_list,list(coef_temp))
  }

  # reshape results
  ## assess_df
  assess_df=data.frame(t(matrix(assess_df,ncol = length(K_seq))));
  rownames(assess_df)=paste("2021",method,K_seq,sep = "-")
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
  return(list(assessment=assess_df,BIC=BIC_seq,eBIC=eBIC_seq,coefficients=coefficient_list,labels=label_df,label_trace=label_trace,fpca=pca_res,reg_dt=reg_dt,
              beta_trace=beta_trace,residuals_trace=residuals_trace,residuals_list=residuals_list,beta_fd=beta_fd_list))
}

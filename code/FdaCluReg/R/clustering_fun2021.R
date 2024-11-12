#' comparison algorithm in 2021
#'
#' @param reg_dt N*M principle components scores of the observations
#' @param label N*1 vector of initialized label sequence
#' @param K number of clusters
#' @param method loss function used in the regression
#'
#' @return list(label result,estimated coefficients,label record matrix)

clustering_fun2021=function(reg_dt,label,K,method="ls")
{
  temp_K_seq=1:K
  label=factor(label,levels = temp_K_seq)
  label_pre=label
  label_after=label
  N=nrow(reg_dt)
  p=ncol(reg_dt)#number of regression variables including constant variable(intercept)
  label_record=rep(1,N)
  beta_trace=NULL
  residuals_trace=NULL

  i=0
  while(i<1|length(which(as.character(label_pre)==as.character(label_after)))<N & i<30)
  {
    label_pre=label_after# previous label_after is label_pre now
    i=i+1#round times
    table_label=table(label_after)
    temp_K_seq=1:K#as.numeric(names(table_label)[table_label>0])
    estimation_temp_round=estimation_fun(reg_dt = reg_dt,label = label_pre,mask_index =NULL,method = method,temp_K_seq=temp_K_seq)
    parameters_temp_round=estimation_temp_round$parameters_mat
    #记录每一轮更新后的beta估计-beta trace
    #empty_parameter_mat=matrix(rep(0,ncol(reg_dt)*K),ncol = K);colnames(empty_parameter_mat)=1:K
    #empty_parameter_mat[,colnames(empty_parameter_mat)%in%colnames(parameters_temp_round)]=parameters_temp_round
    beta_trace=c(beta_trace,list(parameters_temp_round))


    reg_x=as.matrix(cbind(rep(1,N),reg_dt[,-1]))# predictors
    y=reg_dt[,1]# response
    estimated_y=reg_x%*%parameters_temp_round# estimated y for different groups
    temp_residuals=y-estimated_y#
    colnames(temp_residuals)=temp_K_seq
    label_after=colnames(temp_residuals)[apply(abs(temp_residuals),1,which.min)];label_after=factor(label_after,levels = temp_K_seq)
    #记录残差
    residuals_trace=c(residuals_trace,list(temp_residuals))
    #table(label_after)
    #table_label=table(label_after)
    #审查是否观测大于p,如果不是，重新划分类簇
    #abandon_cluster=temp_K_seq[which(table_label<p)]
    #recheck_obs=which(label_after%in%abandon_cluster)
    #if(length(abandon_cluster)>0)
    #{
    #  recheck_abs_residuals=abs(temp_residuals[recheck_obs,!colnames(temp_residuals)%in%abandon_cluster])
    #  if(length(recheck_obs)==1){recheck_abs_residuals=t(as.matrix(recheck_abs_residuals))}
    #  label_after[recheck_obs]=colnames(recheck_abs_residuals)[apply(recheck_abs_residuals,1,which.min)]
    #  }
    #table(label_pre);table(label_after)
    #if(length(unique(label_after))<K){label_after=sample(K,N,replace = T)}#don't need reinitialization
    label_record=cbind(label_record,as.character(label_after))
  }
  #trace 名字
  names(beta_trace)=paste("round",1:length(beta_trace))
  names(residuals_trace)=paste("round",1:length(residuals_trace))

  #回归参数补全
  return(list(label=as.character(label_after),cofficients=estimation_temp_round$parameters_mat,label_record=label_record[,-1],residuals_temp=temp_residuals,beta_trace=beta_trace,residuals_trace=residuals_trace))
}

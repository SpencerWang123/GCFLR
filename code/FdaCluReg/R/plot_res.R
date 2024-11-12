#' plot function for gibbs algorithm and its comparision
#'
#' @param temp_res_gibbs result of gibbs algorithm
#' @param temp_res_2021 result of 2021-comparison algorithm
#' @param K_seq int. vector contains the numbers of clusters
#' @param real_label N*1 label vector contains ground truth labels
#' @param real_coef_mat real coefficient matrix
#' @param y1 ylim for the 1st plot
#' @param y2 ylim for the 2nd plot
#' @param x6 symmetric xlim for the 6th plot
#' @param y6 ylim for the 6th plot
#' @param plot_index logiacal sequence indicating which plot to plot
#'
#' @return automatically generates the plots demanded in the plot_index

plot_res=function(function_list,temp_res_gibbs,temp_res_2021,K_seq,real_label,real_coef_mat,y1=-2,y2=2,y4=c(-10,10),y5=c(-10,10),x6=5,y6=5,plot_index=c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE))
{
  #index1-xi&y_plot;
  #index2-gibbs coef vs real;
  #index3-2021 coef vs real;
  #index4-2021 coef trace;
  #index5-gibbs coef trace;
  #index6-gibbs coef trace density
  color_seq=c("red","blue","orange","purple","green","grey","#cd5c5c","#006400","#FF00FF","#d2691e","#00FFFF","#808000","yellow","black","#FF4000")
  temp_fpca=temp_res_gibbs$fpca[[1]]
  M=max(min(which(cumsum(temp_fpca$varprop)>0.9)))# principle components

  temp_coef_list_gibbs=temp_res_gibbs$coefficients
  temp_coef_list_2021=temp_res_2021$coefficients
  fun1=function_list[[1]]
  fun2=function_list[[2]]
  ##################################plot#####################################
  #visualization of x-y
  #data
  temp_pca=temp_res_gibbs$fpca[[1]]
  pc_scores=temp_pca$scores
  yi_temp=yi_temp
  #plot
  if(plot_index[1]){
    par(mfrow=c(1,M))
    for (m in 1:M)
    {
      label_seq=unique(real_label)
      plot(pc_scores[real_label==label_seq[1],m],yi_temp[real_label==label_seq[1]],col=color_seq[1],main = "xi1-Y",xlab = "X",ylab = "Y")
      for (l in 2:length(label_seq))
      {
        points(pc_scores[real_label==label_seq[l],m],yi_temp[real_label==label_seq[l]],col=color_seq[l])
      }
    }
  }


  ###################################gibbs coef vs real####################################
  #j=1
  if(plot_index[2]){
    par(mfrow=c(1,length(K_seq)))
    for (j in 1:length(K_seq))
    {
      K=K_seq[j]
      temp_type="gibbs"
      ylim_temp=c(-0.5,2)
      temp_coef_list=temp_coef_list_gibbs
      par(new=F)
      temp_coef_mat=temp_coef_list[grep(K,names(temp_coef_list))][[1]]
      pca_fd=temp_fpca$harmonics
      for (k in 1:K)
      {
        beta_estimated=0#intercept
        for (m in 1:M) {beta_estimated=beta_estimated+temp_coef_mat[m+1,k]*pca_fd[m]}
        plot(beta_estimated,xlab = "t",ylab = "y",col=color_seq[k],
             main=paste("beta",temp_type,sep = "-",K,"clusters"),ylim=c(y1,y2),
             xlim = c(0,10));par(new=T)}
      seq_plot=seq(0,10,0.05)
      for (r in 1:nrow(real_coef_mat))
      {
        temp_coef=real_coef_mat[r,]  #one row one cluster
        lines(seq_plot,temp_coef[1]*fun_1(seq_plot)+temp_coef[2]*fun_2(seq_plot),type="l",main=paste("beta",temp_type,sep = "-"),
              ylab = "y",xlab = "t",lty=2);
      }
      #true coef
    }
  }
  ###################################2021 coef vs real####################################
  if(plot_index[3]){
    par(mfrow=c(1,length(K_seq)))
    for (j in 1:length(K_seq))
    {
      K=K_seq[j]
      temp_type="2021"
      ylim_temp=c(-0.5,2)
      temp_coef_list=temp_coef_list_2021
      par(new=F)
      temp_coef_mat=temp_coef_list[grep(K,names(temp_coef_list))][[1]]
      pca_fd=temp_fpca$harmonics
      for (k in 1:K)
      {
        beta_estimated=0#intercept
        for (m in 1:M) {beta_estimated=beta_estimated+temp_coef_mat[m+1,k]*pca_fd[m]}
        plot(beta_estimated,xlab = "t",ylab = "y",col=color_seq[k],
             main=paste("beta",temp_type,sep = "-",K,"clusters"),ylim=c(y1,y2),
             xlim = c(0,10));par(new=T)}
      seq_plot=seq(0,10,0.05)
      for (r in 1:nrow(real_coef_mat))
      {
        temp_coef=real_coef_mat[r,]  #one row one cluster
        lines(seq_plot,temp_coef[1]*fun_1(seq_plot)+temp_coef[2]*fun_2(seq_plot),type="l",main=paste("beta",temp_type,sep = "-"),
              ylab = "y",xlab = "t",lty=2);
      }
      #true coef
    }
  }
  ###################################2021 coef trace####################################
  if(plot_index[4]){
    # comparision beta trace
    # beta trace
    M=M#number of xi
    par(mfrow=c(1,length(K_seq)))
    for (j in 1:length(K_seq))
    {
      K=K_seq[j]
      temp_beta_trace=temp_res_2021$beta_trace[grep(K,names(temp_coef_list))][[1]]
      trace_vec=NULL
      for(l in 1:length(temp_beta_trace))
      {trace_vec=c(trace_vec,as.vector(temp_beta_trace[[l]][-1,]))}
      trace_mat=matrix(trace_vec,ncol = length(temp_beta_trace))

      yu=y4[2];yl=y4[1]
      temp_round=1:ncol(trace_mat)
      for (k in 1:K)
      {
        for (m in 1:(M))
        {
          plot(temp_round,trace_mat[(k-1)*M+m,],col=color_seq[k],ylim = c(yl,yu),
               xlab = "t",ylab = "b",main = paste("2021-K=",K))
          lines(temp_round,trace_mat[(k-1)*M+m,],col=color_seq[k],ylim = c(yl,yu))
          par(new=T)
          par(new=T)
        }}
      par(new=F)
    }

  }

  ###################################gibbs coef trace####################################
  # gibbs beta trace
  # beta trace
  if(plot_index[5]){
    M=M#number of xi
    par(mfrow=c(1,length(K_seq)))
    for (j in 1:length(K_seq))
    {
      K=K_seq[j]
      temp_beta_trace=temp_res_gibbs$beta_trace[grep(K,names(temp_coef_list))][[1]]
      trace_vec=NULL
      for(l in 1:length(temp_beta_trace))
      {trace_vec=c(trace_vec,as.vector(temp_beta_trace[[l]][-1,]))}
      trace_mat=matrix(trace_vec,ncol = length(temp_beta_trace))

      yu=y5[2];yl=y5[1]
      temp_round=1:ncol(trace_mat)
      for (k in 1:K)
      {
        for (m in 1:(M))
        {
          plot(temp_round,trace_mat[(k-1)*M+m,],col=color_seq[k],ylim = c(yl,yu),
               xlab = "t",ylab = "b",main = paste("2021-K=",K))
          lines(temp_round,trace_mat[(k-1)*M+m,],col=color_seq[k],ylim = c(yl,yu))
          par(new=T)
          par(new=T)
        }}
      par(new=F)
    }
  }

  ###################################gibbs coef density####################################
  if(plot_index[6]){
    M=M#number of xi
    par(mfrow=c(1,length(K_seq)))
    for (j in 1:length(K_seq))
    {
      K=K_seq[j]
      temp_beta_trace=temp_res_gibbs$beta_trace[grep(K,names(temp_coef_list))][[1]]
      trace_vec=NULL
      for(l in 1:length(temp_beta_trace))
      {trace_vec=c(trace_vec,as.vector(temp_beta_trace[[l]][-1,]))}
      trace_mat=matrix(trace_vec,ncol = length(temp_beta_trace))

      temp_round=1:ncol(trace_mat)
      for (k in 1:K)
      {
        for (m in 1:(M-1))
        {
          plot(density(trace_mat[(k-1)*M+m,]),col=color_seq[k],xlim=c(-x6,x6),ylim = c(0,y6),
               xlab = "t",ylab = "b",main = paste("gibbs-K=",K))
          par(new=T)
        }}
      par(new=F)
    }}
############BIC###########
  par(mfrow=c(1,1))
  plot(as.numeric(gsub(".+-","",colnames(temp_res_2021$labels))),temp_res_2021$BIC,type = "l",col="blue",xlab = "clusters",ylab = "BIC")
  lines(as.numeric(gsub(".+-","",colnames(temp_res_gibbs$labels))),temp_res_gibbs$BIC,type = "l",col="red")

}

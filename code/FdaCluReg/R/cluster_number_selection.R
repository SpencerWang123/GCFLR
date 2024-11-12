#' cluster_number_selection
#'
#' @param temp_res_gibbs gibbs clustering regression result
#'
#' @return a list consists:
#' height_record : length(K_seq)*max(K_seq) matrix each row records the height of hclust result
#' dist_matrix : list with same length as K_seq whose every element consist p dist matrix indicating the distance between the beta(t) of clusters
#' hclust_trace : hclust results

cluster_number_selection=function(temp_res_gibbs)
  {
  #计算函数型数据系数之间的距离
  fd_list=temp_res_gibbs$beta_fd
  cluster_seq=as.numeric(gsub(".+ ","",names(fd_list)))#cluster_number K_seq
  #变量数
  var_number=ncol(temp_res_gibbs$beta_fd[[1]][[1]]$coefs)

  #空的合并距离表存储和并距离的信息-一会儿的empty_height矩阵
  height_record=vector("list",var_number)
  empty_hight=matrix(rep(0,max(cluster_seq-1)*length(cluster_seq)),nrow = length(cluster_seq))
  for (v in 1:var_number){height_record[[v]]=empty_hight}

  #空的list存储层次聚类结果
  hclust_record_all=vector("list",length = length(cluster_seq))
  names(hclust_record_all)=paste("cluster",cluster_seq)

  #存储不同变量不同类簇的距离矩阵
  cluster_dist_list=vector("list",length = length(cluster_seq))
  names(cluster_dist_list)=paste("cluster num of",cluster_seq)

  for(i in 1:length(fd_list))#逐个K_seq中的类别进行
  {

    temp_fd=fd_list[[i]]
    K=cluster_seq[i]#temp cluster number
    temp_fd_dist=NULL
    for (j in 1:length(temp_fd))
    {
      for(l in 1:length(temp_fd))#fd dist for var_num variables
      {
        temp_gap=temp_fd[[j]]-temp_fd[[l]]# distance fd object
        temp_dist=diag(as.matrix(inprod(temp_gap,temp_gap)))
        temp_fd_dist=c(temp_fd_dist,temp_dist)
        var_num=length(temp_dist)
      }
    }
    temp_fd_dist=matrix(temp_fd_dist,nrow = var_num)#transform the dist into dist matrix
    dist_mat_list=vector("list",var_num)

    #执行聚类
    ##存储聚类结果
    hclust_record=vector("list",var_number)
    names(height_record)=paste("var",1:var_number)

    for (v in 1:var_number)
    {
      temp_dist_mat=matrix(temp_fd_dist[v,],ncol = K)
      dist_mat_list[[v]]=temp_dist_mat
      #hclustering the beta fd
      temp_dist=as.dist(temp_dist_mat)
      hclust_temp=hclust(temp_dist)
      #记录height
      height_record[[v]][i,1:length(hclust_temp$height)]=rev(hclust_temp$height)
      #记录完整结果
      hclust_record[[v]]=hclust_temp
    }
    names(dist_mat_list)=paste("var",1:var_num)
    names(hclust_record)=paste("var",1:var_num)
    #saving the dist mat
    cluster_dist_list[[i]]=dist_mat_list
  }

  return(list(height_record=height_record,dist_matrix=cluster_dist_list,hclust_trace=hclust_record_all))
}


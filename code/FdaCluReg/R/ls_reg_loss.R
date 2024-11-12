#' residuals of the regression using least square error (LS)
#'
#' @param reg_dt N*M matrix indicating the M principle scores of N functional samples
#'
#' @return N*1 vector of residuals of the regression

ls_reg_loss=function(reg_dt,p){
  if(is.null(nrow(reg_dt))){return(rep(0,0))}else
  {
    if(nrow(reg_dt)<p){return(reg_dt[,1]-mean(reg_dt[,1]))}else
    {return(lm(y~.,data=reg_dt)$residuals)}}
}

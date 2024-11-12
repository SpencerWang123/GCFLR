#' residuals of the regression using huber loss function (HL)
#'
#' @param reg_dt N*M matrix indicating the M principle scores of N functional samples
#'
#' @return N*1 vector of residuals of the regression
hl_reg_loss=function(reg_dt,p)
{
  if(is.null(nrow(reg_dt))){return(rep(0,0))}else
  {
    if(nrow(reg_dt)<p){return(reg_dt[,1]-mean(reg_dt[,1]))}else
    {return(MASS::rlm(y~.,data=reg_dt,psi=psi.huber,k=1.345,maxit = 200)$residuals)}}
}

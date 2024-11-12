#' regression using least square error(LSE)
#'
#' @param reg_dt N*M matrix indicating the M principle scores of N functional samples
#'
#' @return (M+1)*1 vector indicating the regression coefficient whose first elements indicates the intercept term

ls_reg_fun=function(reg_dt,p)
{
  if(is.null(nrow(reg_dt))){return(rep(0,p))}else
  {
      if(nrow(reg_dt)<p){return(c(mean(reg_dt[,1]),rep(0,p-1)))}else
      {return(lm(y~.,data=reg_dt)$coefficients)}}
    }




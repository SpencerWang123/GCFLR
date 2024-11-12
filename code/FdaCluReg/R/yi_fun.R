#' generating simulated Y_i
#'
#' @param function_list K function objected aligned in a list
#' @param fda_coef_df_list coefficient matrix with K columns and N rows
#' @param reg_coef_df_list N*K coefficient matrix(true coefficient)
#' @param reg_intercept scalar indicating the intercept of the regression
#' @param reg_var N*1 vector indicating variance parameter of random error term
#' @param reg_error N*1 vector indicating random error
#' @param seed
#'
#' @return N*1 response variable Y with error term

yi_fun=function(function_list,fda_coef_df_list,reg_coef_df_list,reg_intercept,reg_var,reg_error,interval,seed=1)
{
  set.seed(seed)
  N=length(reg_var)# number of observations
  Y_all=reg_intercept+reg_var*reg_error# intercept_term + error_term
  # calculate var by var

  for (var_index in 1:length(fda_coef_df_list))
  {
    #paramters of temp variable
    fda_coef_df=fda_coef_df_list[var_index][[1]]
    reg_coef_df=reg_coef_df_list[var_index][[1]]
    int_term=NULL#save intergral

    for (i in 1:nrow(fda_coef_df))#sample by sample generating
    {
      # x_i(t)
      x=function(t)
      {
        res=0
        for (l in 1:length(function_list))
        {
          temp_function=function_list[l][[1]]
          res=res+fda_coef_df[i,l]*temp_function(t)
        }
        return(res)
      }
      # beta_i(t)
      beta=function(t)
      {
        res=0
        for (l in 1:length(function_list))
        {
          temp_function=function_list[l][[1]]
          res=res+reg_coef_df[i,l]*temp_function(t)
        }
        return(res)
      }
      # intergrate
      fun_int=function(t){return(x(t)*beta(t))}
      intergrate_temp=integrate(fun_int,interval[1],interval[2])$value

      # int_term_i
      int_term=c(int_term,intergrate_temp)
    }
    Y_all=Y_all+int_term# intercept_term + error_term + intergration_term
  }
  return(Y_all)
}

'''
Base code from:

Robert J. DiNapoli, Scott M. Fitzpatrick, Matthew F. Napolitano, Torben C. Rick, Jessica H. Stone, Nicholas P. Jew,
Marine reservoir corrections for the Caribbean demonstrate high intra- and inter-island variability in local reservoir offsets,
Quaternary Geochronology, Volume 61, 2021
https://doi.org/10.1016/j.quageo.2020.101126.

'''

deltar_Can <- read.csv("./deltar_Can.csv", sep=";")

w_mean_fun <- function (df, delta_r, delta_r_error){
  pooled_mean <- sum(delta_r/(delta_r_error^2))/sum(1/(delta_r_error^2)) #error weighted pooled mean
  weighted_uncertainty <- sqrt(1/sum(1/(delta_r_error^2)))#pooled sd
  t_chi_sq <- sum(((delta_r-pooled_mean)^2)/(delta_r_error^2)) #t value
  degree_freedom <- nrow(df)-1
  t_crit <- qchisq(p=0.05, df=degree_freedom, lower.tail=F)
  normalized_chi_sq <- t_chi_sq/(nrow(df)-1)
  est_se <- weighted_uncertainty*sqrt(nrow(df)) #estimated standard error
  ext_var <- sqrt((sd(delta_r)^2)-(est_se^2)) #external variance
  T_uncertainty <- sqrt((weighted_uncertainty^2)+(ext_var^2)) #delta r with external variance
  if(normalized_chi_sq > 1){
    outliers <- "outlier(s)"
  } else {
    outliers <- "no outliers"
  }
  return(data.frame(deparse(substitute(df)), pooled_mean, weighted_uncertainty,
                    t_chi_sq, degree_freedom, t_crit, normalized_chi_sq, outliers, T_uncertainty))
}

w_mean_Can<-w_mean_fun(df= deltar_Can, delta_r = deltar_Can$deltar, delta_r_error = deltar_Can$deltar_error)

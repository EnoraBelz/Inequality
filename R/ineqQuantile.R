#' Means of the ranges according to a GB2 distribution
#'
#' @param p numeric, vector of probabilities
#' @param q numeric, vector of quantiles
#' @param mu numeric, vector of location parameter values
#' @param sigma numeric, vector of scale parameter values
#' @param nu numeric, vector of skewness parameter values
#' @param tau numeric, vector of kurtosis parameter values
#' @seealso \code{\link{compute_LC}}
#' @export
mean_binned_GB2 = function(p,q, mu = 1, sigma = 1, nu = 1, tau = 0.5){

  n=length(q)
  Q1=c(0,q[1:n])
  Q2=c(q[1:n],Inf)

  P1=c(0,p[1:n])
  P2=c(p[1:n],1)
  E=rep(NA,n-1)

  f = function(x) (x * dGB2(x, mu, sigma, nu, tau))/(p2 - p1)


  for(i in 1:n){
    q1=Q1[i]
    q2=Q2[i]
    p1=pGB2(q1,mu,sigma,nu,tau)
    p2=pGB2(q2,mu,sigma,nu,tau)
    E[i]=integrate(f,lower=q1,upper=q2)$value
  }

  q1= Q1[n+1]
  q2 = Inf
  p1= pGB2(q1, mu, sigma, nu, tau)
  p2= pGB2(q2, mu, sigma, nu, tau)
  itg=try(integrate(f,q1,q2)$value,silent=TRUE)
  if(inherits(itg, "try-error")) itg=try(integrate(f,q1,1e8, subdivisions = 100000L)$value,silent=TRUE)
  if(inherits(itg, "try-error")) itg=try(integrate(f,q1,1e7, subdivisions = 100000L)$value,silent=TRUE)
  E[n+1]=itg

  names(E)=paste("(",100*P1,",",100*P2,")",sep="")
  return(E)
}

#' Means of the bins according the conditional expectation or midpoint methods
#'
#' @param ID vector of ID area
#' @param p numeric, vector of probabilities
#' @param bound_min numeric, minimum of the interval
#' @param bound_max numeric, maximum of the interval
#' @param nb numeric, number of the interval
#' @param method string, type of methods ("CondExp" for conditional expectation method of "Midpoint" for midpoint method)
#' @param whichpareto numeric, probability from which pareto tail is assumed
#' @return A dataframe with the bounds, the mean, the cumulative income share and the cumulative population share.
#' @seealso \code{\link{run_compute_LC}}
#' @references Belz (2019), \emph{Estimating Inequality Measure from Quantile Data}
#' @examples
#' data("tabulated_income")
#' BelAir5 = tabulated_income[tabulated_income$iris=="Bel Air 5",]
#' compute_LC(ID=BelAir5$iris, p=BelAir5$prop_cum_population, bound_min = BelAir5$bound_min, bound_max = BelAir5$bound_max, nb = BelAir5$prop_population, method = "CondExp")
#' compute_LC(ID=BelAir5$iris, p=BelAir5$prop_cum_population, bound_min = BelAir5$bound_min, bound_max = BelAir5$bound_max, nb = BelAir5$prop_population, method = "Midpoint")
#' @export
compute_LC = function(ID, p, bound_min, bound_max, nb, method="CondExp", whichpareto=.8){
  n = length(ID)

  if(method=="CondExp"){
    x = run_GB_family(ID=ID,
                      hb = nb,
                      bin_min = bound_min,
                      bin_max = bound_max,
                      obs_mean =rep(NA, length(ID)),
                      ID_name = "State",
                      q = seq(0.006,0.996, length.out = 1000),
                      modelsToFit = c('GB2'))

    P = p[1:(n-1)] #probability vector
    Q = bound_max[1:(n-1)] #quantile vector

    mean <- mean_binned_GB2(p=P,
                            q=Q,
                            mu=as.numeric(x$params$GB2[1]),
                            sigma=as.numeric(x$params$GB2[2]),
                            nu=as.numeric(x$params$GB2[3]),
                            tau=as.numeric(x$params$GB2[4])) }

  if(method=="Midpoint"){
    #fit Pareto
    Q90 = bound_min[n]
    p90 = nb[n]
    A = bound_max[which(p==whichpareto)]
    a = log(p90)/(log(A)-log(Q90))
    mean_open = Q90*(a)/(a-1)

    mean = c((bound_min[1:(n-1)]+bound_max[1:(n-1)])/2,mean_open)
  }

  # Calcul des parts cumules
  total = nb*mean
  income_cum = cumsum(total)/sum(total)
  population_cum = p

  df = data.frame(ID,bound_min,bound_max,mean,income_cum,population_cum)
  return(df)}


#' Means of the bins according the conditional expectation or midpoint methods for several areas
#'
#' @param ID vector of ID area
#' @param p numeric, vector of probabilities
#' @param bound_min numeric, minimum of the interval
#' @param bound_max numeric, maximum of the interval
#' @param nb numeric, number of the interval
#' @param method string, type of methods ("CondExp" for conditional expectation method of "Midpoint" for midpoint method)
#' @param whichpareto numeric, probability from which pareto tail is assumed
#' @return A dataframe with the bounds, the mean, the cumulative income share and the cumulative population share.
#' @references Belz (2019), \emph{Estimating Inequality Measure from Quantile Data}
#' @examples
#' data("tabulated_income")
#' run_compute_LC(ID=tabulated_income$iris, p=tabulated_income$prop_cum_population, bound_min = tabulated_income$bound_min, bound_max = tabulated_income$bound_max, nb = tabulated_income$prop_population, method = "CondExp")
#' @export
run_compute_LC = function(ID,p, bound_min, bound_max, nb, method="CondExp", whichpareto=.8){
  data <- data.frame(ID, p, bound_min, bound_max, nb)
  area = unique(ID)
  res = NA
  for(i in 1:length(area)){
    data_to_fit = data %>% filter(ID == area[i])
    res=rbind(res,
              compute_LC(ID = data_to_fit$ID,
                         p = data_to_fit$p,
                         bound_min = data_to_fit$bound_min,
                         bound_max = data_to_fit$bound_max,
                         nb = data_to_fit$nb,
                         method = method,
                         whichpareto = whichpareto))
  }

  return(res %>% filter(is.na(ID)==F))}


#' Functional form Kakwani and Podder (1973)
#' @param p numeric, vector of probabilities
#' @param a numeric
#' @param b numeric
#' @return The value of the Lorenz curve at the point \code{p}
#' @seealso \code{\link{RGKO}}, \code{\link{Arnold}}, \code{\link{Ortega}}, \code{\link{Chotikapanich}}, \code{\link{Sarabia}}, \code{\link{Rohde}}
#' @references Kakwani and Podder (1973), \emph{On the estimation of Lorenz curves from grouped observations}
#' @examples
#' KP(0.5, 0.9, 1.5)
#' 0.253135
#' @export
KP = function(p,a,b) p^a*exp(-b*(1-p))


#' Functional form  Rasche et al. (1980)
#' @param p numeric, vector of probabilities
#' @param a numeric
#' @param b numeric
#' @return The value of the Lorenz curve at the point \code{p}
#' @seealso \code{\link{KP}}, \code{\link{Arnold}}, \code{\link{Ortega}}, \code{\link{Chotikapanich}}, \code{\link{Sarabia}}, \code{\link{Rohde}}
#' @references  Rasche et al. (1980), \emph{Functional forms for estimating the Lorenz curve: comment}
#' @examples
#' RGKO(0.5, 0.7, 1.5)
#' 0.2383538
#' @export
RGKO = function(p,a,b) (1-(1-p)^a)^b


#' Functional form  Arnold (1986)
#' @param p numeric, vector of probabilities
#' @param a numeric
#' @param b numeric
#' @return The value of the Lorenz curve at the point \code{p}
#' @seealso \code{\link{KP}}, \code{\link{RGKO}}, \code{\link{Ortega}}, \code{\link{Chotikapanich}}, \code{\link{Sarabia}}, \code{\link{Rohde}}
#' @references  Arnold (1986), \emph{A class of hyperbolic Lorenz curves}
#' @examples
#' Arnold(0.5, 1.2, 2.2)
#' 0.3333333
#' @export
Arnold = function(p,a,b) (p*(1 + (a-1)*p))/(1 + (a-1)*p + b*(1-p))


#' Functional form  Ortega et al. (1991)
#' @param p numeric, vector of probabilities
#' @param a numeric
#' @param b numeric
#' @return The value of the Lorenz curve at the point \code{p}
#' @seealso \code{\link{KP}}, \code{\link{RGKO}}, \code{\link{Arnold}},  \code{\link{Chotikapanich}}, \code{\link{Sarabia}}, \code{\link{Rohde}}
#' @references  Ortega et al. (1991), \emph{A new functional form for estimating Lorenz curves}
#' @examples
#' Ortega(0.5, 0.5, 0.7)
#' 0.2718315
#' @export
Ortega = function(p,a,b) p^a*(1 - (1-p)^b)

#' Functional form  Chotikapanich (1993)
#' @param p numeric, vector of probabilities
#' @param k numeric
#' @return The value of the Lorenz curve at the point \code{p}
#' @seealso \code{\link{KP}}, \code{\link{RGKO}}, \code{\link{Arnold}},  \code{\link{Ortega}}, \code{\link{Sarabia}}, \code{\link{Rohde}}
#' @references  Chotikapanich (1993), \emph{A comparison of alternative functional forms for the Lorenz curve}
#' @examples
#' Chotikapanich(0.5, 2.2)
#' 0.2497399
#' @export
Chotikapanich = function(p,k) (exp(k*p)-1)/(exp(k)-1)


#' Functional form  Sarabia (1997)
#' @param p numeric, vector of probabilities
#' @param pi1 numeric
#' @param pi2 numeric
#' @param a1 numeric
#' @param a2 numeric
#' @return The value of the Lorenz curve at the point \code{p}
#' @seealso \code{\link{KP}}, \code{\link{RGKO}}, \code{\link{Arnold}},  \code{\link{Ortega}}, \code{\link{Chotikapanich}}, \code{\link{Rohde}}
#' @references  Sarabia (1997), \emph{A hierarchy of Lorenz curves based on the generalized Tukey’s lambda distribution}
#' @examples
#' Sarabia(0.5, 0.1, 0.7, 1.8, 0.3)
#' 0.2885717
#' @export
Sarabia = function(p, pi1, pi2, a1, a2) pi1*p + pi2*p^a1 + (1 - pi1 - pi2)*(1 - (1-p)^a2)

#' Functional form  Rohde (2009)
#' @param p numeric, vector of probabilities
#' @param b numeric
#' @return The value of the Lorenz curve at the point \code{p}
#' @seealso \code{\link{KP}}, \code{\link{RGKO}}, \code{\link{Arnold}},  \code{\link{Ortega}}, \code{\link{Chotikapanich}}, \code{\link{Sarabia}}
#' @references  Rohde (2009), \emph{An alternative functional form for estimating the Lorenz curve}
#' @examples
#' Rohde(0.5, 1.5)
#' 0.25
#' @export
Rohde = function(p, b){p*((b-1)/(b-p))}


#' Optimisation of a parametric Lorenz curve
#' @param ID vector of ID area
#' @param income_cum numeric, vector of cumulaive income shares
#' @param population_cum numeric, vector of cumulative population shares
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @return A dataframe with functional form, the parameters, the value of the NLS and the Chi-squared statistic
#' @seealso \code{\link{run_optim_LC}}
#' @export
optim_LC = function(ID,income_cum, population_cum, function_form){

  # likelihood
  logLik = function(theta){
    F=function(x){LOI(x,theta)}
    LCtheo = F(population_cum)
    NLS = sum((LCtheo-income_cum)^2)
    return(list(NLS=NLS))
  }

  if(function_form == "KP"){
    LOI = function(x,theta) KP(x, a = theta[1], b = theta[2])
    # constraints
    ui = matrix(c(1,0,0,1),2,2)
    ci = c(1,0)
    theta0 = c(1.5,1.25)
  }
  if(function_form == "RGKO"){
    LOI = function(x,theta) RGKO(x, a = theta[1], b = theta[2])
    # constraints
    ui = matrix(c(1,0,-1,0,0,1),3,2,byrow= T)
    ci = c(0,-1,1)
    theta0 = c(0.5,2)
  }
  if(function_form == "ARNOLD"){
    LOI = function(x,theta) Arnold(x, a = theta[1], b = theta[2])
    # constraints
    ui = matrix(c(1,0,-1,+1,0,1),3,2,byrow=T)
    ci = c(0,-1,0)
    theta0 = c(1,1)
  }
  if(function_form == "CHOTIKAPANICH"){
    LOI = function(x,theta) Chotikapanich(x, k = theta[1])
    # constraints
    ui = matrix(c(1),1,1)
    ci = c(0)
    theta0 = c(2)
  }
  if(function_form == "SARABIA"){
    LOI = function(x,theta) Sarabia(x, pi1= theta[1], pi2= theta[2], a1= theta[3], a2 = theta[4])
    # constraints
    ui = matrix(c(1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1,0,0,0,-1),5,4,byrow=T)
    ci = c(0,-1,1,0,-1)
    theta0 = c(0.5,0.5,1.5,0.5)
  }
  if(function_form == "ORTEGA"){
    LOI = function(x,theta) Ortega(x, a = theta[1], b = theta[2])
    # constraints
    ui = matrix(c(1,0,0,1,0,-1),3,2,byrow= T)
    ci = c(0,0,-1)
    theta0 = c(1,0.3)
  }
  if(function_form == "ROHDE"){
    LOI = function(x,theta) Rohde(x, b = theta[1])
    # constraints
    ui = matrix(c(1),1,1)
    ci = c(1)
    theta0 = c(2)
  }

  opt_chisq = constrOptim(theta=theta0, f=function(x) logLik(x)$NLS, grad=NULL,  ui=ui, ci=ci)
  par = opt_chisq$par
  NLS = opt_chisq$value
  LCtheo = LOI(population_cum,theta=par)
  chisq = sum((income_cum-LCtheo)^2/LCtheo)

  return(data.frame(ID = unique(ID),function_form = function_form, par1 = par[1],par2 = par[2],par3 = par[3],par4 = par[4],
                    NLS=NLS,chisq=chisq))
}



#' Optimisation of a parametric Lorenz curve for several functional forms
#' @param ID vector of ID area
#' @param income_cum numeric, vector of cumulaive income shares
#' @param population_cum numeric, vector of cumulative population shares
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @return A dataframe with functional form, the parameters, the value of the NLS and the Chi-squared statistic
#' @seealso \code{\link{run_optim_LC}}
#' @export
optim_LC_function = function(ID,income_cum, population_cum, function_form){
  zone = unique(ID)
  res = as.data.frame(t(sapply(function_form,
                               function(i) optim_LC(ID = ID,
                                                    income_cum = income_cum ,
                                                    population_cum = population_cum,
                                                    function_form = i))))
  res = res %>% mutate(ID = zone,
                       function_form = rownames(res))
  return(res)
}


#' Optimisation of a parametric Lorenz curve for several functional forms and several areas
#' @param ID vector of ID area
#' @param income_cum numeric, vector of cumulaive income shares
#' @param population_cum numeric, vector of cumulative population shares
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @return A dataframe with functional form, the parameters, the value of the NLS and the Chi-squared statistic
#' @examples
#' data("tabulated_income")
#' LC_tabulated_income = run_compute_LC(ID=tabulated_income$ID, p=tabulated_income$prop_cum_population, bound_min = tabulated_income$bound_min, bound_max = tabulated_income$bound_max, nb = tabulated_income$prop_population, method = "CondExp")
#' run_optim_LC(ID = unique(LC_tabulated_income$ID),income_cum = LC_tabulated_income$income_cum, population_cum=LC_tabulated_income$population_cum, function_form = c("KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA", "ROHDE"))
#' @export
run_optim_LC = function(ID,income_cum, population_cum, function_form){
  data <- data.frame(ID,income_cum, population_cum)
  area = unique(ID)
  res = NA
  for(i in 1:length(area)){
    data_to_fit = data %>% filter(ID == area[i])
    res=rbind(res,
              optim_LC_function(ID = data_to_fit$ID,
                                income_cum = data_to_fit$income_cum,
                                population_cum  = data_to_fit$population_cum,
                                function_form = function_form))
    print(paste(area[i]))
  }
  res = res %>% filter(is.na(ID)==F)
  res = res %>% mutate(par1=as.numeric(par1),
                                           par2=as.numeric(par2),
                                           par3=as.numeric(par3),
                                           par4=as.numeric(par4),
                                           NLS = as.numeric(NLS),
                                           chisq = as.numeric(chisq))

  return(res)}


#' Gini of the functional form  Kakwani and Podder (1973)
#' @param a numeric
#' @param b numeric
#' @return The value of the Gini index
#' @seealso \code{\link{Gini_RGKO}}, \code{\link{Gini_Arnold}}, \code{\link{Gini_Ortega}},  \code{\link{Gini_Chotikapanich}}, \code{\link{Gini_Sarabia}}, \code{\link{Gini_Rohde}}
#' @references  Kakwani and Podder (1973), \emph{On the estimation of Lorenz curves from grouped observations}
#' @examples
#' Gini_KP(0.9, 1.5)
#' 0.3331993
#' @export
Gini_KP = function(a,b) as.numeric(1 - 2*exp(-b)*fAsianOptions::kummerM(b,1+a,2+a)/(1+a))


#' Gini of the functional form  Rasche et al. (1980)
#' @param a numeric
#' @param b numeric
#' @return The value of the Gini index
#' @seealso \code{\link{Gini_KP}}, \code{\link{Gini_Arnold}}, \code{\link{Gini_Ortega}},  \code{\link{Gini_Chotikapanich}}, \code{\link{Gini_Sarabia}}, \code{\link{Gini_Rohde}}
#' @references  Rasche et al. (1980), \emph{Functional forms for estimating the Lorenz curve: comment}
#' @examples
#' Gini_RGKO(0.7, 1.5)
#' 0.3868911
#' @export
Gini_RGKO = function(a,b) 1 - 2/a * beta(1/a,b + 1)


#' Gini of the functional form  Arnold (1986)
#' @param a numeric
#' @param b numeric
#' @return The value of the Gini index
#' @seealso \code{\link{Gini_KP}}, \code{\link{Gini_RGKO}}, \code{\link{Gini_Ortega}},  \code{\link{Gini_Chotikapanich}}, \code{\link{Gini_Sarabia}}, \code{\link{Gini_Rohde}}
#' @references  Arnold (1986), \emph{A class of hyperbolic Lorenz curves}
#' @examples
#' Gini_Arnold(1.2, 2.2)
#' 0.3484886
#' @export
Gini_Arnold = function(a,b) b/(b-a+1) + 2*a*b/(b-a+1)^2*(1 + (b+1)/(b-a+1)*log(a/(b+1)))

#' Gini of the functional form  Ortega et al. (1991)
#' @param a numeric
#' @param b numeric
#' @return The value of the Gini index
#' @seealso \code{\link{Gini_KP}}, \code{\link{Gini_RGKO}}, \code{\link{Gini_Arnold}},  \code{\link{Gini_Chotikapanich}}, \code{\link{Gini_Sarabia}}, \code{\link{Gini_Rohde}}
#' @references  Ortega et al. (1991),  \emph{A new functional form for estimating Lorenz curves}
#' @examples
#' Gini_Ortega(0.5, 0.7)
#' 0.3310822
#' @export
Gini_Ortega = function(a,b) (a-1)/(a+1) + 2*beta(a+1,b+1)


#' Gini of the functional form  Chotikapanich (1993)
#' @param k numeric
#' @return The value of the Gini index
#' @seealso \code{\link{Gini_KP}}, \code{\link{Gini_RGKO}}, \code{\link{Gini_Arnold}},  \code{\link{Gini_Ortega}}, \code{\link{Gini_Sarabia}}, \code{\link{Gini_Rohde}}
#' @references  Chotikapanich (1993),  \emph{A comparison of alternative functional forms for the Lorenz curve}
#' @examples
#' Gini_Chotikapanich(2.2)
#' 0.3401299
#' @export
Gini_Chotikapanich = function(k) ((k-2)*exp(k) + (k+2))/(k*(exp(k)-1))

#' Gini of the functional form  Sarabia (1997)
#' @param pi1 numeric
#' @param pi2 numeric
#' @param a1 numeric
#' @param a2 numeric
#' @return The value of the Gini index
#' @seealso \code{\link{Gini_KP}}, \code{\link{Gini_RGKO}}, \code{\link{Gini_Arnold}},  \code{\link{Gini_Ortega}}, \code{\link{Gini_Chotikapanich}}, \code{\link{Gini_Rohde}}
#' @references  Sarabia (1997), \emph{A hierarchy of Lorenz curves based on the generalized Tukey’s lambda distribution}
#' @examples
#' Gini_Sarabia(0.1, 0.7, 1.8, 0.3)
#' 0.3076923
#' @export
Gini_Sarabia = function(pi1, pi2, a1, a2 ) pi2*(1-2/(1+a1)) - (1-pi1-pi2)*(1-2/(1+a2))


#' Gini of the functional form  Rohde (2009)
#' @param p numeric, vector of probabilities
#' @param b numeric
#' @return The value of the Gini index
#' @seealso \code{\link{Gini_KP}}, \code{\link{Gini_RGKO}}, \code{\link{Gini_Arnold}},  \code{\link{Gini_Ortega}}, \code{\link{Gini_Chotikapanich}}, \code{\link{Gini_Sarabia}}
#' @references  Rohde (2009), \emph{An alternative functional form for estimating the Lorenz curve}
#' @examples
#' Gini_Rohde(1.5)
#' 0.3520816
#' @export
Gini_Rohde = function(b) 2*b*((b-1)*log((b-1)/b)+1) - 1



#' Computation of Gini index according to a functional form and its parameters
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @param par1 numeric, parameter of the functional form
#' @param par2 numeric, parameter of the functional form
#' @param par3 numeric, parameter of the functional form
#' @param par4 numeric, parameter of the functional form
#' @return The value of the Gini index accorting a functional form and its parameters
#' @seealso \code{\link{compute_Pietra}}, \code{\link{compute_TH}}, \code{\link{compute_TL}}, \code{\link{compute_topshare}}
#' @examples
#' data("tabulated_income")
#' LC_tabulated_income = run_compute_LC(ID=tabulated_income$ID, p=tabulated_income$prop_cum_population, bound_min = tabulated_income$bound_min, bound_max = tabulated_income$bound_max, nb = tabulated_income$prop_population, method = "CondExp")
#' Optim_LC = run_optim_LC(ID = unique(LC_tabulated_income$ID),income_cum = LC_tabulated_income$income_cum, population_cum=LC_tabulated_income$population_cum, function_form = c("KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA", "ROHDE"))
#' Optim_LC
#' compute_Gini(function_form = "KP", par1 = 1, par2 = 1.450156)
#' 0.3488301
#' Optim_LC %>% ungroup() %>%  rowwise() %>% mutate(Gini = compute_Gini(function_form, par1, par2, par3, par4)) %>% select(ID, function_form, Gini)
#' @export
compute_Gini = function(function_form, par1,par2 = NA,par3 = NA,par4 = NA){
  if (function_form=="RGKO"){Gini_est = Gini_RGKO(par1,par2)}
  if (function_form=="ARNOLD"){Gini_est = Gini_Arnold(par1,par2)}
  if (function_form=="CHOTIKAPANICH"){Gini_est=Gini_Chotikapanich(par1)}
  if (function_form=="ORTEGA"){Gini_est = Gini_Ortega(par1,par2)}
  if (function_form=="SARABIA"){Gini_est = Gini_Sarabia(par1,par2,par3,par4)}
  if (function_form=="ROHDE"){Gini_est = Gini_Rohde(par1)}
  if (function_form=="KP"){Gini_est = Gini_KP(par1,par2)}
  return(Gini_est)
}


#' Computation of Pietra index according to a functional form and its parameters
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @param par1 numeric, parameter of the functional form
#' @param par2 numeric, parameter of the functional form
#' @param par3 numeric, parameter of the functional form
#' @param par4 numeric, parameter of the functional form
#' @return The value of the Pietra index accorting a functional form and its parameters
#' @seealso \code{\link{compute_Gini}}, \code{\link{compute_TH}}, \code{\link{compute_TL}}, \code{\link{compute_topshare}}
#' @examples
#' data("tabulated_income")
#' LC_tabulated_income = run_compute_LC(ID=tabulated_income$ID, p=tabulated_income$prop_cum_population, bound_min = tabulated_income$bound_min, bound_max = tabulated_income$bound_max, nb = tabulated_income$prop_population, method = "CondExp")
#' Optim_LC = run_optim_LC(ID = unique(LC_tabulated_income$ID),income_cum = LC_tabulated_income$income_cum, population_cum=LC_tabulated_income$population_cum, function_form = c("KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA", "ROHDE"))
#' Optim_LC
#' compute_Pietra(function_form = "KP", par1 = 1, par2 = 1.450156)
#' 0.2645622
#' Optim_LC %>% ungroup() %>%  rowwise() %>% mutate(Pietra = compute_Pietra(function_form, par1, par2, par3, par4)) %>% select(ID, function_form, Pietra)
#' @export
compute_Pietra = function(function_form, par1, par2, par3=NA, par4=NA){
  p = seq(0,1,0.000001)
  if (function_form=="RGKO"){L_p = RGKO(p,par1,par2)}
  if (function_form=="ARNOLD"){L_p  = Arnold(p,par1,par2)}
  if (function_form=="CHOTIKAPANICH"){L_p = Chotikapanich(p,par1)}
  if (function_form=="ORTEGA"){L_p  = Ortega(p,par1,par2)}
  if (function_form=="SARABIA"){L_p = Sarabia(p,par1,par2,par3,par4)}
  if (function_form=="ROHDE"){L_p = Rohde(p,par1)}
  if (function_form=="KP"){L_p = KP(p,par1,par2)}
  Pietra_est = max(p - L_p)
  return(Pietra_est)
}


# Lprime:
# Compute the derivative of a function L at the point x
# @x : point
# @L : function
# @h
#' Computation of the derivative of a function L at the point x
#' @param x numeric, point
#' @param L function
#' @param h numeric
#' @export
Lprime=function(x,L,h=1e-5){
  d= sapply(x, function(x)
    if(x < h)  (L(x+h)-L(x))/h
    else if (x>(1-h)) (L(x)-L(x-h))/h
    else if ((x>=h)&(x<=(1-h))) (L(x+h)-L(x-h))/(2*h))
  return(d)}


#' Computation of Theil'L index according to a functional form and its parameters
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @param par1 numeric, parameter of the functional form
#' @param par2 numeric, parameter of the functional form
#' @param par3 numeric, parameter of the functional form
#' @param par4 numeric, parameter of the functional form
#' @return The value of the Theil'L index accorting a functional form and its parameters
#' @seealso \code{\link{compute_Gini}}, \code{\link{compute_Pietra}}, \code{\link{compute_TL}}, \code{\link{compute_topshare}}
#' @examples
#' data("tabulated_income")
#' LC_tabulated_income = run_compute_LC(ID=tabulated_income$ID, p=tabulated_income$prop_cum_population, bound_min = tabulated_income$bound_min, bound_max = tabulated_income$bound_max, nb = tabulated_income$prop_population, method = "CondExp")
#' Optim_LC = run_optim_LC(ID = unique(LC_tabulated_income$ID),income_cum = LC_tabulated_income$income_cum, population_cum=LC_tabulated_income$population_cum, function_form = c("KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA", "ROHDE"))
#' Optim_LC
#' compute_TL(function_form = "KP", par1 = 1, par2 = 1.450156)
#' 0.2109571
#' Optim_LC %>% ungroup() %>%  rowwise() %>% mutate(TL = compute_TL(function_form, par1, par2, par3, par4)) %>% select(ID, function_form, TL)
#' @export
compute_TL = function(function_form, par1, par2, par3=NA, par4=NA){
  if (function_form=="RGKO"){L_p = function(p) RGKO(p,par1,par2)}
  if (function_form=="ARNOLD"){L_p  = function(p) Arnold(p,par1,par2)}
  if (function_form=="CHOTIKAPANICH"){L_p = function(p) Chotikapanich(p,par1)}
  if (function_form=="ORTEGA"){L_p  = function(p) Ortega(p,par1,par2)}
  if (function_form=="SARABIA"){L_p = function(p) Sarabia(p,par1,par2,par3,par4)}
  if (function_form=="ROHDE"){L_p = function(p) Rohde(p,par1)}
  if (function_form=="KP"){L_p = function(p) KP(p,par1,par2)}

  TL = integrate(f=function(x) -log(Lprime(x,L=L_p,h=1e-5)),lower=0,upper=1)

  return(TL$value)
}

#' Computation of Theil'H index according to a functional form and its parameters
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @param par1 numeric, parameter of the functional form
#' @param par2 numeric, parameter of the functional form
#' @param par3 numeric, parameter of the functional form
#' @param par4 numeric, parameter of the functional form
#' @return The value of the Theil'H index accorting a functional form and its parameters
#' @seealso \code{\link{compute_Gini}}, \code{\link{compute_Pietra}}, \code{\link{compute_TL}}, \code{\link{compute_topshare}}
#' @examples
#' data("tabulated_income")
#' LC_tabulated_income = run_compute_LC(ID=tabulated_income$ID, p=tabulated_income$prop_cum_population, bound_min = tabulated_income$bound_min, bound_max = tabulated_income$bound_max, nb = tabulated_income$prop_population, method = "CondExp")
#' Optim_LC = run_optim_LC(ID = unique(LC_tabulated_income$ID),income_cum = LC_tabulated_income$income_cum, population_cum=LC_tabulated_income$population_cum, function_form = c("KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA", "ROHDE"))
#' Optim_LC
#' compute_TH(function_form = "KP", par1 = 1, par2 = 1.450156)
#' 0.1900282
#' Optim_LC %>% ungroup() %>%  rowwise() %>% mutate(TH = compute_TH(function_form, par1, par2, par3, par4)) %>% select(ID, function_form, TH)
#' @export
compute_TH = function(function_form, par1, par2, par3=NA, par4=NA){
  if (function_form=="RGKO"){L_p = function(p) RGKO(p,par1,par2)}
  if (function_form=="ARNOLD"){L_p  = function(p) Arnold(p,par1,par2)}
  if (function_form=="CHOTIKAPANICH"){L_p = function(p) Chotikapanich(p,par1)}
  if (function_form=="ORTEGA"){L_p  = function(p) Ortega(p,par1,par2)}
  if (function_form=="SARABIA"){L_p = function(p) Sarabia(p,par1,par2,par3,par4)}
  if (function_form=="ROHDE"){L_p = function(p) Rohde(p,par1)}
  if (function_form=="KP"){L_p = function(p) KP(p,par1,par2)}

  TH = integrate(f=function(x) Lprime(x,L=L_p,h=1e-5)*log(Lprime(x,L=L_p,h=1e-5)),lower=0,upper=1)

  return(TH$value)
}

#' Computation of topshare according to a functional form and its parameters
#' @param p numeric, value of the topshare
#' @param function_form string, functional form in "KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA" and "ROHDE"
#' @param par1 numeric, parameter of the functional form
#' @param par2 numeric, parameter of the functional form
#' @param par3 numeric, parameter of the functional form
#' @param par4 numeric, parameter of the functional form
#' @return The value of the topshare accorting a functional form and its parameters
#' @seealso \code{\link{compute_Gini}}, \code{\link{compute_Pietra}}, \code{\link{compute_TL}}, \code{\link{compute_TH}}
#' @examples
#' data("tabulated_income")
#' LC_tabulated_income = run_compute_LC(ID=tabulated_income$ID, p=tabulated_income$prop_cum_population, bound_min = tabulated_income$bound_min, bound_max = tabulated_income$bound_max, nb = tabulated_income$prop_population, method = "CondExp")
#' Optim_LC = run_optim_LC(ID = unique(LC_tabulated_income$ID),income_cum = LC_tabulated_income$income_cum, population_cum=LC_tabulated_income$population_cum, function_form = c("KP", "RGKO", "ARNOLD", "CHOTIKAPANICH", "SARABIA", "ORTEGA", "ROHDE"))
#' Optim_LC
#' compute_topshare(p = 0.95, function_form = "KP", par1 = 1, par2 = 1.450156)
#' 0.1164444
#' Optim_LC %>% ungroup() %>%  rowwise() %>% mutate(topshare95= compute_topshare(p=0.95, function_form, par1, par2, par3, par4)) %>% select(ID, function_form, topshare95)
#' @export
compute_topshare = function(p, function_form, par1, par2, par3=NA, par4=NA){
  if (function_form=="RGKO"){L_p = RGKO(p,par1,par2)}
  if (function_form=="ARNOLD"){L_p  = Arnold(p,par1,par2)}
  if (function_form=="CHOTIKAPANICH"){L_p = Chotikapanich(p,par1)}
  if (function_form=="ORTEGA"){L_p  = Ortega(p,par1,par2)}
  if (function_form=="SARABIA"){L_p = Sarabia(p,par1,par2,par3,par4)}
  if (function_form=="ROHDE"){L_p = Rohde(p,par1)}
  if (function_form=="KP"){L_p = KP(p,par1,par2)}
  return(topshare = 1 - L_p)
}


#' import HMisc
#' import knitr
#' importFrom("tidyverse", "binequality")

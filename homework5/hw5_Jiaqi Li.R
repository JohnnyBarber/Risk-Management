library(pbivnorm)
library(dplyr)
#------------------------------------3---------------------------------------#
Greeks = function(S1,S2,K,r,sigma1,sigma2,tau,rho){
  gamma1 = (log(S1/K)+(r-0.5*sigma1^2)*tau)/(sigma1*sqrt(tau))
  gamma2 = (log(S2/K)+(r-0.5*sigma2^2)*tau)/(sigma2*sqrt(tau))
  sigma = sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2)
  a1 = gamma1+sigma1*sqrt(tau)
  b1 = (log(S2/S1)-0.5*sigma^2*tau)/(sigma*sqrt(tau))
  theta1 = (rho*sigma2-sigma1)/sigma
  a2 = gamma2+sigma2*sqrt(tau)
  b2 = (log(S1/S2)-0.5*sigma^2*tau)/(sigma*sqrt(tau))
  theta2 = (rho*sigma1-sigma2)/sigma
  out = as.data.frame(matrix(c(a1,a2,b1,b2,theta1,theta2,gamma1,gamma2),
                             nrow = 2, ncol = 4, byrow = F))
  names(out) = c("alpha","beta","theta","gamma")
  return(out)
}


OptionPrice = function(r,S1_0,S2_0,sigma1,sigma2,K,rho,tau,mu1,mu2){
  
  ABTG = Greeks(S1_0,S2_0,K,r,sigma1,sigma2,tau,rho)
  
  N2_sig1 = pbivnorm(ABTG[1,1],ABTG[1,2],rho = ABTG[1,3])
  N2_sig2 = pbivnorm(ABTG[2,1],ABTG[2,2],rho = ABTG[2,3])
  N2_gam = pbivnorm(ABTG[1,4],ABTG[2,4],rho = 0.4)
  
  M = S1_0*N2_sig1+S2_0*N2_sig2-K*exp(-r*tau)*N2_gam
  return(M)
}

r = 0.00005*252; S1_0 = 99; S2_0 = 101; sigma1 = sigma2 = 0.02*sqrt(252)
K = 100; rho = 0.4; tau = 0.5; mu1 = mu2 = 0.0003*252

M = OptionPrice(r,S1_0,S2_0,sigma1,sigma2,K,rho,tau,mu1,mu2)
M

#------------------------------------4---------------------------------------#
dS1 = 99.01; dS2 = 101.01; dS3 = 99-0.01; dS4 = 101-0.01

M_dS1 = OptionPrice(r,dS1,S2_0,sigma1,sigma2,K,rho,tau,mu1,mu2)
M_dS2 = OptionPrice(r,S1_0,dS2,sigma1,sigma2,K,rho,tau,mu1,mu2)
M_dS3 = OptionPrice(r,dS3,S2_0,sigma1,sigma2,K,rho,tau,mu1,mu2)
M_dS4 = OptionPrice(r,S1_0,dS4,sigma1,sigma2,K,rho,tau,mu1,mu2)
M_pp = OptionPrice(r,dS1,dS2,sigma1,sigma2,K,rho,tau,mu1,mu2)

delta1 = (M_dS1-M)/0.01
delta2 = (M_dS2-M)/0.01

gamma1 = (M_dS1-2*M+M_dS3)/0.01^2
gamma2 = (M_dS2-2*M+M_dS4)/0.01^2
gamma12 = (M_pp-M_dS1-M_dS2+M)/0.01^2

#delta approach
VaR_delta = (-(delta1*mu1/252*S1_0+delta2*mu2/252*S2_0)+
  qnorm(0.99,0,1)*sqrt((delta1*sigma1/sqrt(252)*S1_0)^2+
                         (delta2*sigma2/sqrt(252)*S2_0)^2+
                         2*delta1*delta2*rho*sigma1*sigma2/252*S1_0*S2_0))/M

#delta-gamma approach
VaR_delta_gamma = (-(delta1*mu1/252*S1_0+delta2*mu2/252*S2_0+
          0.5*gamma1*sigma1^2/252*S1_0^2+0.5*gamma2*sigma2^2/252*S2_0^2+
          gamma12*sigma1*sigma2/252*S1_0*S2_0*rho)+
  qnorm(0.99,0,1)*sqrt((delta1*sigma1/sqrt(252)*S1_0)^2+(delta2*sigma2/sqrt(252)*S2_0)^2+
                     2*delta1*delta2*sigma1*sigma2/252*S1_0*S2_0*rho))/M

#------------------------------------5---------------------------------------#
StockPrice = function(S1_0,S2_0,r,sd1,sd2,mu1,mu2,rho,Time,paths,steps){
  dt = Time/steps
  Z = rnorm(paths*steps*2,0,1)
  Z1 = Z[1:(steps*paths)] %>% matrix(nrow = paths, ncol = steps) %>% as.data.frame()
  Z2 = Z[(paths*steps+1):(2*steps*paths)] %>% matrix(nrow = paths, ncol = steps) %>% as.data.frame()
  dW1t = sd1*Z1
  dW2t = sd2*rho*Z1+sd2*sqrt(1-rho^2)*Z2
  S1t = matrix(0,nrow = paths, ncol = steps+1)
  S2t = matrix(0,nrow = paths, ncol = steps+1)
  S1t[,1] = S1_0
  S2t[,1] = S2_0
  for(i in 1:steps){
    S1t[,i+1] = S1t[,i]*exp((mu1-0.5*sd1^2)*dt+dW1t[,i]*sqrt(dt))
    S2t[,i+1] = S2t[,i]*exp((mu2-0.5*sd2^2)*dt+dW2t[,i]*sqrt(dt))
  }
  return(list(S1 = S1t, S2 = S2t))
}

r = 0.00005*252; S1_0 = 99; S2_0 = 101; sigma1 = sigma2 = 0.02*sqrt(252)
K = 100; rho = 0.4; tau = 0.5; mu1 = mu2 = 0.0003*252
steps = as.integer(252/2); paths = 1000

set.seed(1234)

test1 = StockPrice(S1_0,S2_0,r,sigma1,sigma2,mu1,mu2,rho,tau,paths,steps)[[1]]
test2 = StockPrice(S1_0,S2_0,r,sigma1,sigma2,mu1,mu2,rho,tau,paths,steps)[[2]]
Price0 = cbind(test1[,steps+1],test2[,steps+1])

Min_stock0 = apply(Price0,1,min)
temp0 = cbind(Min_stock0-K,0)
option_Mt0 = apply(temp0,1,max)*exp(-r*tau)
Mt0 = mean(option_Mt0)

Mt1 = c()
for(i in 1:1000){
  test3 = StockPrice(test1[,2][i],test2[,2][i],r,sigma1,sigma2,mu1,mu2,rho,tau-tau/steps,paths,steps-1)[[1]]
  test4 = StockPrice(test1[,2][i],test2[,2][i],r,sigma1,sigma2,mu1,mu2,rho,tau-tau/steps,paths,steps-1)[[2]]
  Price1 = cbind(test3[,steps],test4[,steps])
  Min_stock1 = apply(Price1,1,min)
  temp1 = cbind(Min_stock1-K,0)
  option_Mt1 = apply(temp1,1,max)*exp(-r*(tau-tau/steps))
  Mt1[i] = mean(option_Mt1)
}

VaR_sim = -quantile((Mt1-Mt0)/Mt0,0.01)

VaR_delta
VaR_delta_gamma
VaR_sim

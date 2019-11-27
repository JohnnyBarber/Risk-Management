library(dplyr)
library(ggplot2)
library(data.table)

####################################################################
#                           Question 1                             #
####################################################################
data = read.csv("hw3_return.csv") %>% as.data.table()
data = data[,Date := as.Date(as.character(Date),format = "%m/%d/%Y")]
data = data[, Year := as.numeric(format(Date, "%Y"))]
setkey(data,Year)
test = data[.(2015:2017)]

n = length(data[.(2014)]$Return)
N = length(data$Return)
c = 0.99
f_VaR = function(c,n,N,data){
  VaR = c()
  for(i in n:N){
    VaR[i-n+1] = -(quantile(data[1:i],1-c,na.rm = T))
  }
  return(VaR)
}
VaR = f_VaR(c,n,N,data$Return)
test[, VaR := VaR]
qplot(Date,Return,data = test, geom = "line", color = I("steelblue")) + 
  theme_bw() + geom_line(aes(y=-VaR), color = I("red"))
exceptions = length(data[Return<(-VaR),]$Return)

lambda = 0.995
f_exp_VaR = function(N,c,lambda,data){
  weights = c()
  for(i in 1:N){
    weights[i] = lambda^(N-i)*(1-lambda)/(1-lambda^(N))
  }
  test = as.data.table(cbind(Return = data$Return[1:N],weights))
  setorder(test,Return)
  temp = cumsum(test$weights)
  test[,cum_sum := temp]
  VaR_exp = test[cum_sum > 1-c,]$Return[1]
  return(VaR_exp)
}

VaR_exp = c()
for(i in n:N){
  VaR_exp[i-n+1] = f_exp_VaR(i,c,lambda,data)
}
test[,VaR_exp := VaR_exp]

qplot(Date,Return,data = test, geom = "line", color = I("steelblue")) + 
  geom_line(aes(y=VaR_exp), color = I("red")) + theme_bw()
exceptions = length(test[Return<VaR_exp,]$Return)

####################################################################
#                           Question 2                             #
####################################################################

mu = c();std = c()
for(i in 1:length(test$Return)){
  mu[i] = mean(test$Return[1:i])
  std[i] = sd(test$Return[1:i])
}
x = mu-2.33*std
fx = dnorm(x,mean = mu, sd = std)
StdDev = (1/fx)*sqrt((0.99*0.01)/n)
CI_left = -test$VaR-abs(qnorm(0.025,0,1))*StdDev
CI_right = -test$VaR+abs(qnorm(0.025,0,1))*StdDev
qplot(Date,Return,data = test, geom = "line", color = I("steelblue")) + 
  theme_bw() + geom_line(aes(y=-test$VaR), color = I("red"))+
  geom_line(aes(y=CI_left), color = I("red"),linetype = 2) +
  geom_line(aes(y=CI_right), color = I("red"),linetype = 2) +
  geom_ribbon(aes(ymin=CI_left,ymax=CI_right), fill="red", alpha="0.2") 

bootstrap_VaR_H = c()
bootstrap_VaR_L = c()
for(i in 1:(N-n)){
  bootstrap_VaR = c()
  for(j in 1:1000){
    temp_bootstrap = sample(data$Return[1:(n+i-1)],
                                size = n+i-1, replace = T)
    bootstrap_VaR[j] = quantile(temp_bootstrap,0.01,type = 1)
  }
  bootstrap_VaR_H[i] = quantile(bootstrap_VaR,0.975)
  bootstrap_VaR_L[i] = quantile(bootstrap_VaR,0.025)
}

qplot(Date,Return,data = test, geom = "line", color = I("steelblue")) + 
  theme_bw() + geom_line(aes(y=-VaR), color = I("red"))+
  geom_line(aes(y=bootstrap_VaR_L), color = I("red"),linetype = 2) +
  geom_line(aes(y=bootstrap_VaR_H), color = I("red"),linetype = 2) +
  geom_ribbon(aes(ymin=bootstrap_VaR_L,ymax=bootstrap_VaR_H), fill="red", alpha="0.2") 


bootstrap_expVaR_H = c()
bootstrap_expVaR_L = c()
for(i in 1:(N-n)){
  bootstrap_expVaR = c()
  for(j in 1:500){
    temp_bootstrap_exp = sample(data$Return[1:(n+i-1)],
                            size = n+i-1, replace = T)
    weights = lambda^(n+i-1-c(1:(n+i-1)))*(1-lambda)/(1-lambda^(n+i-1))
    temp_exp = cbind(temp_bootstrap_exp,weights) %>% as.data.table()
    setorder(temp_exp,temp_bootstrap_exp)
    x = cumsum(temp_exp$weights)
    temp_exp[,cum_sum := x]
    bootstrap_expVaR[j] = temp_exp[cum_sum>1-c,]$temp_bootstrap_exp[1]
  }
  bootstrap_expVaR_H[i] = quantile(bootstrap_expVaR,0.975)
  bootstrap_expVaR_L[i] = quantile(bootstrap_expVaR,0.025)
}

qplot(Date,Return,data = test, geom = "line", color = I("steelblue")) + 
  theme_bw() + geom_line(aes(y=VaR_exp), color = I("red"))+
  geom_line(aes(y=bootstrap_expVaR_L), color = I("red"),linetype = 2) +
  geom_line(aes(y=bootstrap_expVaR_H), color = I("red"),linetype = 2) +
  geom_ribbon(aes(ymin=bootstrap_expVaR_L,ymax=bootstrap_expVaR_H),
              fill="red", alpha="0.2") 


####################################################################
#                           Question 3                             #
####################################################################
data[, `:=` (Year = as.numeric(format(Date,"%Y")), Month = as.numeric(format(Date,"%m")))]
data[,volatility := sd(Return), by = list(Year,Month)]
data[,volatility_shift := shift(volatility)]
data[,norm_RET := Return/volatility_shift]
hist(data$norm_RET,breaks = 30,col=rgb(0,1,0,1/4))
hist(data$Return,breaks = 30,col=rgb(0,0,1,1/4))


####################################################################
#                           Question 4                             #
####################################################################
VaR_new = f_VaR(c,n,N,data$norm_RET)[-1]
setkey(data,Year)
test = na.omit(data[.(2015:2017)])
test[, VaR_new := VaR_new]
qplot(Date,norm_RET,data = test, geom = "line", color = I("steelblue")) + 
  theme_bw() + geom_line(aes(y=-VaR_new), color = I("red"))
exceptions = length(test[norm_RET<(-VaR_new),]$Return)
exceptions

test = data[.(2014:2015)]
f_exp_VaR = function(N,c,lambda,data){
  weights = c()
  for(i in 1:N){
    weights[i] = lambda^(N-i)*(1-lambda)/(1-lambda^(N))
  }
  test = as.data.table(cbind(Return = data$norm_RET[1:N],weights))
  setorder(test,Return)
  temp = cumsum(test$weights)
  test[,cum_sum := temp]
  VaR_exp = test[cum_sum > 1-c,]$Return[1]
  return(VaR_exp)
}
VaR_exp = c()
for(i in n:N){
  VaR_exp[i-n+1] = f_exp_VaR(i,c,lambda,data)
}
test[,VaR_exp := VaR_exp]

qplot(Date,norm_RET,data = test, geom = "line", color = I("steelblue")) + 
  geom_line(aes(y=test$VaR_exp), color = I("red")) + theme_bw()
exceptions = length(test[norm_RET<VaR_exp,]$Return)
exceptions

####################################################################
#                           Question 5                             #
####################################################################

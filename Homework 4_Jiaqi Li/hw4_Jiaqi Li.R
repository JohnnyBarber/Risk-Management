library(data.table)
library(ggplot2)
library(dplyr)
options(digits = 8)

####################################################################
#                           Question 1                             #
####################################################################
data = read.csv("hw3_return.csv") %>% as.data.table()
data[,Date := as.Date(as.character(Date),format="%m/%d/%Y")]
data[,lag_ret := shift(Return)]
data = na.omit(data)
data[,volatilty := sd(lag_ret)]

data[,`:=`(Year = as.integer(format(Date,"%Y")), Month = as.integer(format(Date,"%m")))]
setkey(data,Year)
lambda = 0.94
sigma0 = sd(data[.(2014)]$Return)
Sample = data[.(2015:2017)]
volatility = c(sigma0)
for(i in 1:length(Sample$Return)){
  volatility[i+1] = sqrt(lambda*volatility[i]^2+(1-lambda)*Sample$Return[i]^2)
}
Sample[,estimate := volatility[-1]]
Sample[,VaR := -(qnorm(0.01,0,1)*estimate)]

# curve = as.data.frame(cbind(as.character(rep(Sample$Date,3)),
#                             c(Sample$Return,Sample$estimate,Sample$VaR))) %>%
#   as.data.table()
# curve[,V1 := as.Date(as.character(V1),format="%Y-%m-%d")]
# n = length(Sample$Return)
# curve[,type := c(rep("Return",n),rep("volatility",n),rep("VaR",n))]
# colnames(curve) = c("Date","value","type")

qplot(Date,Return,data = Sample, geom = "line", color = I("steelblue"),
      main = "Return(blue) volatility(yellow) VaR(red)") +
  theme_bw() + geom_line(aes(y=-VaR),color = "firebrick") +
  geom_line(aes(y=estimate),color="darkorange")


####################################################################
#                           Question 2                             #
####################################################################
garch=function(w,a,b,eps){
  l=length(eps)
  sigma_2=c(w/(1-a-b))
  for (i in 1:l){
    sigma_2[i+1]=w+a*eps[i]^2+b*sigma_2[i]
  }
  return (sigma_2)
}

loglike=function(vP,eps=data[.(2015:2017)]$Return){
  w=vP[1]; a=vP[2]; b=vP[3]
  sigma_2=garch(w,a,b,eps)
  logL=-sum(-log(sigma_2)-(eps^2)/sigma_2)
  return (logL)
}

vP0=c(0.2,0.5,0.4)
coeff = optim(fn=loglike,par=vP0)$par

sigma = sqrt(garch(coeff[1],coeff[2],coeff[3],data[.(2015:2017)]$Return))
GARCH_VaR = -(qnorm(0.01,0,1)*sigma)
qplot(Date,Return,data = Sample, geom = "line",color = I("steelblue"))+
  geom_line(aes(y=-GARCH_VaR[-1]),color = "firebrick")+
  geom_line(aes(y=sigma[-1]),color = "darkorange")+theme_bw()


####################################################################
#                           Question 3                             #
####################################################################
Sample[,vol := sd(Return),by=list(Year,Month)]
Sample[,norm_ret := Return/vol]
hist(Sample$norm_ret)

setorder(Sample,norm_ret)
Sample[,loss := -norm_ret]
Sample[,Rank := .I]
Sample[,`:=` (Prob = Rank/length(Sample$Date))]
test = Sample[loss > 0,]
lm.out = lm(log(Prob)~log(loss), data = test)
psy = -1/coef(lm.out)[2]
plot(log(test$loss),log(test$Prob),type="o",xlim=c(-1.0,1.3))

loglike_Pareto = function(para){
  beta = para[1]
  psy = para[2]
  loss = -Sample$norm_ret
  u = quantile(loss,0.95)
  v = loss[loss>u]
  logsum = -sum(log((1+psy*(v-u)/beta)^(-1/psy-1)/beta))
  return(logsum)
}

coeff = optim(fn = loglike_Pareto, par = c(0.8,0.9))$par

loss = -Sample$norm_ret
u = quantile(loss,0.95)
a = quantile(loss,0.99)
N = length(loss)
Nu = sum(loss>=u)
VaR_Pareto = (((N/Nu)*0.01)^(-coeff[2])-1)*coeff[1]/coeff[2]+u


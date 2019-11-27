
####################################################################
#                           Question 2                             #
####################################################################
library(quantmod)
getSymbols(Symbols = "C",from = "2006-01-01",
           to = "2008-12-31", src = "yahoo")

#-------------------------------part a-----------------------------#
library(ggplot2)
library(data.table)
C = as.data.table(C)

compute_VaR = function(C){
  n = length(C)
  PRC = as.data.frame(C)[,n]
  C = cbind(C,PRC)
  C[,Year := format(index,"%Y")]
  setkey(C,Year)
  
  C[, price_shift := shift(PRC)]
  C[, RET := PRC/price_shift-1]
  C = C[-1,]
  C[, Mean := sum(C$RET[1:.I])/.I]
  quant = c()
  for(i in 1:length(C$RET)){
    quant[i] = quantile(C$RET[1:i],0.01)
  }
  C = cbind(C,quant)
  C[, VaR := -(Mean + quant)]
  return(C)
}

C = compute_VaR(C)

qplot(C[.("2008")]$index,C[.("2008")]$VaR, geom = "line",
      xlab = "date", ylab = "VaR", main = "Plot of VaR") + theme_bw()

#-------------------------------part b-----------------------------#
qplot(C[.("2008")]$index,C[.("2008")]$RET, geom = "line",
      xlab = "date", ylab = "VaR", main = "Plot of Backtest") +
  theme_bw() + geom_line(aes(y=-C[.("2008")]$VaR),color = "red")

C[, YesNoMaybe := ifelse(RET < -VaR, 1, 0)]
Exceed = sum(C[.("2008")]$YesNoMaybe)
setkey(C,Year)
n = length(C[.("2008")]$index)
Prob_exceed = Exceed/n

#-------------------------------part c-----------------------------#


####################################################################
#                           Question 3                             #
####################################################################
Bank = c("GS","UBS","JPM","BCS","DB","BAC","BNP.PA","CS","MS","C")
for(i in Bank){
  getSymbols(Symbols = i, from = "2006-01-01",
             to = "2008-12-31", src = "yahoo")
}

#-------------------------------part a-----------------------------#
for(i in 1:length(Bank)){
  assign(Bank[i],compute_VaR(as.data.table(get(Bank[i]))))
}

MS = MS[,list(date = index,RET_MS = RET)]
CS = CS[,list(date = index,RET_CS = RET)]
DB = DB[,list(date = index,RET_DB = RET)]
GS = GS[,list(date = index,RET_GS = RET)]
JPM = JPM[,list(date = index,RET_JPM = RET)]
UBS = UBS[,list(date = index,RET_UBS = RET)]
BAC = BAC[,list(date = index,RET_BAC = RET)]
BCS = BCS[,list(date = index,RET_BCS = RET)]
BNP.PA = BNP.PA[,list(date = index,RET_BNP.PA = RET)]
C = C[,list(date = index,RET_C = RET)]

for(i in 1:length(Bank)){
  setkey(get(Bank[i]),date)
}
Portfolio = merge(GS,UBS,all = T)
for(i in 3:length(Bank)){
  Portfolio = merge(Portfolio,get(Bank[i]),all = T)
}
Portfolio = na.omit(Portfolio)

returns = Portfolio[,-1]

weights = as.matrix(rep(c(1,2),5))
VaR = sapply(returns,function(x){quantile(x,0.01,na.rm = T)})
VaR_p = -sum(VaR*weights)

returns = as.data.table(returns)
weight_M = diag(rep(0.1,10))
DVaR = (-(apply((rep(c(1,2),5)+weight_M)*VaR,2,sum))-VaR_p)/0.1
CVaR = DVaR*weights

profits = apply(returns,2,mean)*weights
RAROC = as.data.frame(profits/CVaR)
rownames(RAROC) = names(returns)
RAROC

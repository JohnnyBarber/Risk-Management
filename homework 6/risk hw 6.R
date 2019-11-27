library(knitr)
library(data.table)

####################################################################
#                           Question 2                             #
####################################################################
#------------------------------part1-------------------------------#
R = 0.6
CDS3 = 0.005
lamb1 = CDS3/(1-R)

CDS5 = 0.006
cal_lamb2 = function(x,lamb1,CDS5,R){
  left = (exp(-3*lamb1)*(1-exp(-2*x))/x+(1-exp(-3*lamb1))/lamb1)*CDS5/(1-R)
  right = 1-exp(-2*x-3*lamb1)
  return(left-right)
}

lamb2 = uniroot(cal_lamb2, lamb1=lamb1,CDS5=CDS5,R=R,interval = c(-1,1))$root

CDS10 = 0.01
cal_lamb3 = function(x,lamb1,lamb2,CDS10,R){
  left = ((1-exp(-3*lamb1))/lamb1+exp(-3*lamb1+3*lamb2)*(exp(-3*lamb2)-exp(-5*lamb2))/lamb2+exp(-3*lamb1-2*lamb2+5*x)*(exp(-5*x)-exp(-10*x))/x)*CDS10/(1-R)
  right = 1 - exp(-3*lamb1-2*lamb2-5*x)
  return(left-right)
}

lamb3 = uniroot(cal_lamb3,lamb1 = lamb1, lamb2 = lamb2, CDS10 = CDS10, R = R, interval = c(-1,1))$root

lambda = c(lamb1,lamb2,lamb3)

kable(t(lambda),col.names=c("0-3","3-5","5-10"),caption="Lambdas")

#------------------------------part2-------------------------------#
c = 0.03*100
survive = c()
cond_default = c()
#no default at all
for(i in 1:(2*6)){
  if(i <= 3*2){
    survive[i] = exp(-i*0.5*lamb1)
  }else if(i > 3*2 & i <= 5*2){
    survive[i] = exp(-3*lamb1-0.5*(i-6)*lamb2)
  }else if(i > 5*2 & i <= 6*2){
    survive[i] = exp(-3*lamb1 - 2*lamb2 - 0.5*(i-10)*lamb3)
  }
}

PV = data.table(cf = c(rep(3,11),103))
PV[, time := (.I)*0.5]
PV = cbind(PV,survive)
PV[, default := 1-survive]
PV[, cond_default := default-shift(default)]
PV[1,"cond_default"] = PV[1,"default"]
PV[, default_pv := R*100+((.I)-1)*c]
PV[, pv := cond_default*default_pv]
PV_bond = sum(PV[,"pv"])+sum(PV[,"cf"])*PV[12,"survive"]

####################################################################
#                           Question 3                             #
####################################################################
#------------------------------part1-------------------------------#
#P(0)
trans_0<-diag(8)
rownames(trans_0)<-c("AAA","AA","A","BBB","BB","B","CCC","Default")
colnames(trans_0)<-c("AAA","AA","A","BBB","BB","B","CCC","Default")
trans_0

#P(1)
trans_1<-matrix(0,nrow=8,ncol=8)
trans_1[1,]=c(90.81, 8.33, 0.68, 0.06, 0.12, 0, 0, 0)
trans_1[2,]=c(0.7, 90.65, 7.79, 0.64, 0.06, 0.14, 0.02, 0)
trans_1[3,]=c(0.09, 2.27, 91.05, 5.52, 0.74, 0.26, 0.01, 0.06)
trans_1[4,]=c(0.02, 0.33, 5.95, 86.93, 5.3, 1.17, 1.12, 0.18)
trans_1[5,]=c(0.03, 0.14, 0.67, 7.73, 80.53, 8.84, 1, 1.06)
trans_1[6,]=c(0, 0.11, 0.24, 0.43, 6.48, 83.46, 4.07, 5.2)
trans_1[7,]=c(0.22, 0, 0.22, 1.3, 2.38, 11.24, 64.86, 19.79)
trans_1[8,]=c(0,0,0,0,0,0,0,100)
rownames(trans_1)<-c("AAA","AA","A","BBB","BB","B","CCC","Default")
colnames(trans_1)<-c("AAA","AA","A","BBB","BB","B","CCC","Default")
trans_1

#------------------------------part6-------------------------------#
moodys = matrix(c(0.000,0.022,0.062,0.174,1.110,3.904,15.894,
                  0.013,0.068,0.199,0.504,3.071,9.274,27.003,
                  0.013,0.136,0.434,0.906,5.371,14.723,35.800,
                  0.037,0.260,0.679,1.373,7.839,19.509,42.796,
                  0.104,0.410,0.958,1.862,10.065,23.869,48.828,
                  0.241,0.682,1.615,2.872,13.911,31.774,56.878,
                  0.489,1.017,2.759,4.632,19.323,40.560,66.212),nrow = 7, ncol = 7)

#------------------------------part7-------------------------------#
M = data.table(coupon = c(rep(3,11),103))
M[,time := 0.5*(.I)]
Q = default["BBB",]
lambda = -(log(1-Q)-log(shift(1-Q)))
lambda[1] = -log(1-Q[1])
lambda = lambda/c(1,1,1,1,1,2,3)

survive = c()
#no default at all
for(i in 1:(2*6)){
  if(i <= 1*2){
    survive[i] = exp(-i*0.5*lambda[1])
  }else if(i > 1*2 & i <= 2*2){
    survive[i] = exp(-lambda[1]-0.5*(i-2)*lambda[2])
  }else if(i > 2*2 & i <= 3*2){
    survive[i] = exp(-lambda[1] - lambda[2] - 0.5*(i-4)*lambda[3])
  }else if(i > 3*2 & i <= 4*2){
    survive[i] = exp(-lambda[1] - lambda[2] - lambda[3] - 0.5*(i-6)*lambda[4])
  }else if(i > 4*2 & i <= 5*2){
    survive[i] = exp(-lambda[1] - lambda[2] - lambda[3] - lambda[4] - 0.5*(i-8)*lambda[5])
  }else if(i > 5*2){
    survive[i] = exp(-lambda[1] - lambda[2] - lambda[3] - lambda[4] - lambda[5] - 0.5*(i-10)*lambda[6])
  }
}
M[, survive := survive]
M[, default := 1-survive]
M[, cond_default := default - shift(default)]
M[1,"cond_default"] = M[1,"default"]
M[, default_pv := R*100+((.I)-1)*c]
M[, pv := cond_default*default_pv]
bond_price = sum(M[,"pv"])+sum(M[,"coupon"])*M[12,"survive"]



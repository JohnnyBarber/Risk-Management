# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 18:14:41 2019

@author: Jiaqi Li
"""
import numpy as np
import scipy.stats as si
import matplotlib.pyplot as plt

def StockPrices(S0, r, sd, T, paths, steps):
    dt = T/steps
    # Generate stochastic process and its antithetic paths
    Z = np.random.normal(0, 1, paths//2 * (steps)).reshape(paths//2,(steps))
    Z_inv = -Z
    dWt = np.sqrt(dt) * Z
    dWt_inv = np.sqrt(dt) * Z_inv
    # bind the normal and antithetic Wt
    dWt = np.concatenate((dWt, dWt_inv), axis=0)
    # define the initial value of St
    St = np.zeros((paths, steps + 1))
    St[:, 0] = S0
    for i in range (steps):
        St[:, i+1] = St[:, i]*np.exp((r - 1/2*(sd**2))*dt + sd*dWt[:, i])
    return St

def f_BS(S0,K,T,sigma,r,type):
    d_1=(np.log(S0/K)+(r+0.5*sigma**2)*T)/(sigma * np.sqrt(T))
    d_2=d_1-sigma*np.sqrt(T)
    if type=="call":
        option=(S0*si.norm.cdf(d_1,0.0,1.0)-K*np.exp(-r*T)*si.norm.cdf(d_2,0.0,1.0))
    if type=="put":
        option=(K*np.exp(-r*T)*si.norm.cdf(-d_2, 0.0, 1.0)-S0*si.norm.cdf(-d_1,0.0,1.0))
    return option
    

###############################################################################
#                                 2.Question 2                                #
###############################################################################
VaR_my = 0.386
VaR_partner = 0.1996
probs = np.arange(0,1.1,0.1)
x = VaR_my*(1-probs)+VaR_partner*probs
plt.rcParams.update({'font.size': 16})
plt.figure(1,figsize=(7,5))
ax = plt.plot(probs,x)
plt.title("VaR with different probabilities on 2 VaRs")
plt.xlabel("probabilities")
plt.ylabel("VaR")

###############################################################################
#                                 3.Question 3                                #
###############################################################################
S0 = 50
r = 0.02
sd = 0.16
T = 1/4 #3 months maturity
paths = 1000 #1000 simulations
steps = 21*3 #91 days in 3 months
K = 50 #ATM call option stick price

stock = StockPrices(S0,r,sd,T,paths,steps) #simulate 10000 stock prices
S10 = stock[:,10] #take stock prices at day 10
C0 = f_BS(S0,K,T,sd,r,"call")
#compute Call optoin price today using black-scholes
C10 = f_BS(S10,K,T-10/252,sd,r,"call")
#compute call option prices at day 10 for all simulations
Cc = np.quantile(C10,0.01)
#take the 0.01th option price at day 10 from the distribution
VaR_C = C0 - Cc #compute VaR

N_Coption = 100/VaR_C #number of options to buy
Cost_C = C0*N_Coption #total cost of all the options bought
Bond_C = 100-Cost_C #total bond to buy/short
w_option_C = Cost_C/(Cost_C+Bond_C) #weight on call option
w_bond_C = 1 - w_option_C #weight on bond


###############################################################################
#                               3.Question 4                                  #
###############################################################################
S0 = 50
r = 0.02
sd = 0.16
T = 1/4 #3 months maturity
paths = 1000 #1000 simulations
steps = 21*3 #91 days in 3 months
K = 50 #ATM put option stick price

stock = StockPrices(S0,r,sd,T,paths,steps) #simulate 10000 stock prices
S10 = stock[:,10] #take stock prices at day 10
P0 = f_BS(S0,K,T,sd,r,"put")
#compute Put optoin price today using black-scholes
P10 = f_BS(S10,K,T-10/252,sd,r,"put")
#compute Put option prices at day 10 for all simulations
Pc = np.quantile(P10,0.99)
#take the 0.99th option price at day 10 from the distribution
VaR_P = P0 - Pc #compute VaR

N_Poption = 100/VaR_P #number of options to buy
Cost_P = P0*N_Poption #total cost of all the options bought
Bond_P = 100-Cost_P #total bond to buy/short
w_option_P = Cost_P/(Cost_P+Bond_P) #weight on put option
w_bond_P = 1 - w_option_P #weight on bond


###############################################################################
#                                 3.Question 5                                #
###############################################################################
S0 = 50
r = 0.02
sd = 0.16
T = 1/4 #3 months maturity
paths = 1000 #1000 simulations
steps = 21*3 #91 days in 3 months
K = np.arange(30,71) #put option stick price range
n = len(K) #length of all put option stick prices

#call option-------------------------------------------------------------------
C0 = np.zeros(n)
for i in range(n):
    C0[i] = f_BS(S0,K[i],T,sd,r,"call")
#compute Put optoin price today using black-scholes
stock = StockPrices(S0,r,sd,T,paths,steps) #simulate 10000 stock prices
S10 = stock[:,10] #take stock prices at day 10
C10 = np.zeros((len(S10),n))
Cc = np.zeros(n)
for i in range(n):
    C10[:,i] = f_BS(S10,K[i],T-10/252,sd,r,"call")
    #compute Put option price at day 10 for each simulation
    Cc[i] = np.quantile(C10[:,i],0.01)
    #take the 0.99th option price at day 10 from the distribution
VaR_C = C0 - Cc #compute VaR

N_Coption = 100/VaR_C #number of options to buy
Cost_C = C0*N_Coption #total cost of all the options bought
Bond_C = 100-Cost_C #total bond to buy/short
w_option_C = Cost_C/(Cost_C+Bond_C) #weight on put option
w_bond_C = 1 - w_option_C #weight on bond


#put option--------------------------------------------------------------------
P0 = np.zeros(n)
for i in range(n):
    P0[i] = f_BS(S0,K[i],T,sd,r,"put")
#compute Put optoin price today using black-scholes
P10 = np.zeros((len(S10),n))
Pc = np.zeros(n)
for i in range(n):
    P10[:,i] = f_BS(S10,K[i],T-10/252,sd,r,"put")
    #compute Put option price at day 10 for each simulation
    Pc[i] = np.quantile(P10[:,i],0.99)
    #take the 0.99th option price at day 10 from the distribution
VaR_P = P0 - Pc #compute VaR

N_Poption = 100/VaR_P #number of options to buy
Cost_P = P0*N_Poption #total cost of all the options bought
Bond_P = 100-Cost_P #total bond to buy/short
w_option_P = Cost_P/(Cost_P+Bond_P) #weight on put option
w_bond_P = 1 - w_option_P #weight on bond

#plot the weights--------------------------------------------------------------
plt.figure(2,figsize=(7,5))
plt.rcParams.update({'font.size': 16})
ax1 = plt.plot(K,w_option_C)
ax2 = plt.plot(K,w_option_P)
plt.legend(["call option","put option"])
plt.xlabel("strik price")
plt.ylabel("absolute value of long/short position")
plt.title("Put/Call Option Position on Different Strike Prices")

#return------------------------------------------------------------------------
C_R = np.zeros(n);P_R = np.zeros(n)
for i in range(n):
    C_R[i] = w_option_C[i]*(np.mean(C10[:,i])-C0[i])/C0[i]+w_bond_C[i]*(r)*10/252
    P_R[i] = w_option_P[i]*(np.mean(P10[:,i])-P0[i])/P0[i]+w_bond_P[i]*(r)*10/252
plt.figure(3,figsize=(7,5))
plt.rcParams.update({'font.size': 16})
ax1 = plt.plot(K,C_R)
ax2 = plt.plot(K,P_R)
plt.legend(["call portfolio","put portfolio"])
plt.xlabel("strik price")
plt.ylabel("portfolio return")
plt.title("Put/Call Portfolio Return on Different Strike Prices")

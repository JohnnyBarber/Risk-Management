import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import norm
import scipy.optimize as optimize

df = pd.read_csv('/Users/leonard/Desktop/409/week3/hw3_returns2.csv')
start_date = df[df['Date'] == '1/2/2015'].index[0]

### Problem 1 Revisiting problem 1 from homework 3
## Q1
lambda_EWMA = 0.94
m = 252 # number of past data
df['Variance_EWMA'] = 0
wt = pd.Series((1-lambda_EWMA)*lambda_EWMA**np.arange(m-1,-1,-1))
df.loc[start_date,'Variance_EWMA'] = np.average(df['Return'][(start_date-m):start_date]**2, weights=wt)
for i in range(start_date+1,len(df['Variance_EWMA'])):
    df.loc[i, 'Variance_EWMA'] = (1-lambda_EWMA) * df.loc[i-1,'Return']**2 + lambda_EWMA * df.loc[i-1,'Variance_EWMA']
df['sigma_EWMA'] = np.sqrt(df['Variance_EWMA'])
df['VaR_EWMA'] = -norm.ppf(0.01) * df['sigma_EWMA']
plt.plot(df[start_date:]['Return'], label='Realized Return')
plt.plot(-df[start_date:]['VaR_EWMA'], label='VaR_EWMA')
plt.legend()
plt.show()
exception_EWMA = sum(df['Return'][start_date:]<-df['VaR_EWMA'][start_date:])

## Q2

def GARCH_Log_likelihood(parameters, data):
    alpha, beta, omega = parameters
    data['Variance_GARCH'] = omega/(1-alpha-beta)
    for i in range(1, len(data['Variance_EWMA'])):
        data.loc[i,'Variance_GARCH'] = omega + alpha*data.loc[i-1,'Return']**2 + beta*data.loc[i-1,'Variance_GARCH']
    Log_likelihood = - sum(-np.log(data['Variance_GARCH']) - data['Return']**2/data['Variance_GARCH'])
    return Log_likelihood

Initial_para = [0.3, 0.3, 0.1]
df['Variance_GARCH'] = 0

# data = df.loc[(start_date-m):(start_date-1)].reset_index(drop=True)
# Parameters = optimize.fmin(GARCH_Log_likelihood, Initial_para, args=(data,),maxiter=1000)
# alpha, beta, omega = Parameters
# df['Variance_GARCH'] = omega/(1-alpha-beta)
# for i in range(1,len(df['Variance_GARCH'])):
#     df.loc[i, 'Variance_GARCH'] = omega + alpha*df.loc[i-1, 'Return']**2 + beta*df.loc[i-1, 'Variance_GARCH']

for i in range(start_date,len(df['Variance_GARCH'])):
    data = df.loc[(i-m):(i-1)].reset_index(drop=True)
    if ((i-2)%50==0):
        Parameters = optimize.fmin(GARCH_Log_likelihood, Initial_para, args=(data,),maxiter=500)
        alpha, beta, omega = Parameters
    data['Variance_GARCH'] = omega/(1-alpha-beta)
    for j in range(1, len(data['Variance_EWMA'])):
        data.loc[j,'Variance_GARCH'] = omega + alpha*data.loc[j-1,'Return']**2 + beta*data.loc[j-1,'Variance_GARCH']
    df.loc[i,'Variance_GARCH'] = omega + alpha*data.loc[len(data)-1,'Return']**2 + beta*data.loc[len(data)-1,'Variance_GARCH']

df['sigma_GARCH'] = np.sqrt(df['Variance_GARCH'])
df['VaR_GARCH'] = -norm.ppf(0.01)*df['sigma_GARCH']
plt.plot(df[start_date:]['Return'], label='Realized Return')
plt.plot(-df[start_date:]['VaR_GARCH'], label='VaR_GARCH')
plt.legend()
plt.show()
exception_GARCH = sum(df['Return'][start_date:]<-df['VaR_GARCH'][start_date:])

## Q3

# Assume invest 1 million
df['Historical_Mean'] = df['Return'].rolling(22).mean()
df['Historical_Mean'] = df['Historical_Mean'].shift(1)
df['Historical_Std'] = df['Return'].rolling(22).std()
df['Historical_Std'] = df['Historical_Std'].shift(1)
df['Normalized_Ret'] = (df['Return'] - df['Historical_Mean'])/df['Historical_Std']
df['Normalized_loss'] = -df['Normalized_Ret']*1000000
df['Rank_loss'] = df['Normalized_loss'].rank(ascending=False)
plt.plot(np.log(df[start_date:]['Normalized_loss']), np.log(df[start_date:]['Rank_loss']), 'ro')
plt.xlabel('log_loss')
plt.ylabel('log_rank')
plt.show()
# It is like a straight line when log_loss is greater than 14, so
# it the left tail of these normalized returns follows a power law

def Pareto_likelihood(parameters,data,u):
    beta, epsilo = parameters
    log_likelihood = -sum(np.log(((1/beta)*(1+(epsilo*(data-u))/beta)**(-1/epsilo - 1))))
    return log_likelihood

Initial_para = [0, 0.05]
df['Loss'] = -df['Return']*1000000
u = np.percentile(df['Loss'].dropna(),95)
# data = df[22:(len(df)-2)].loc[df['Normalized_loss'] > u,'Normalized_loss']
data = df.loc[df['Loss'] > u,'Loss']

beta, epsilo = optimize.fmin(Pareto_likelihood,Initial_para,args=(data, u))
Nu = len(data)
# n = len(df[22:(len(df)-2)])
n = len(df)-1
c = 0.01
VaR_Pareto = u + beta/epsilo*((n/Nu*(1-c))**(-epsilo)-1)

### Problem 2 LTCM
## Q2








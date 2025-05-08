#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 14:25:53 2018

@author: david
"""

import math
import random
from matplotlib.pyplot import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

odds=[2,2]
#proba=[0.4, 0.2, 0.4]
#proba=[0.2, 0.4, 0.4]
proba=[0.2, 0.8]
#bets=[0.5, 0.25, 0.25]

p=proba[0]
o=odds[0]
x=np.sqrt(p*(1-p))
alpha_c=x/(1.0/o-p+x)
alpha=0.4
b=p+(1-alpha)*x/alpha

bets=[b,1-b]

#bets=proba # Kelly
stat=100 # number of races
S=[0] # initial capital
logcapital=0
list=[]

def Kelly(proba,odds):
      return np.sum(np.multiply(proba,np.log(np.multiply(odds,proba))))

for i in range(stat):
         print(i)
         sum_proba=proba[0]
         j=0
         still=1
         r=random.random()
         while still==1: 
             if r<sum_proba:           
                 horse=j
                 still=0
             else:
                 j=j+1
                 sum_proba=sum_proba+proba[j]
         #print('winner',horse+1) 
         list.append(horse)  
#         capital=capital+np.log(odds[horse]*bets[horse])
         logcapital=np.log(odds[horse]*bets[horse])
         S.append(logcapital)

#print(S)
print('Critical alpha',alpha_c)
print('Mean growth rate',np.mean(S))
print('Variance',np.var(S))
print('Kellys prediction',Kelly(proba,odds))
#lambda_estimated=np.log(np.mean(exp))/(T*dt)        
plot(S)
legend()



# Histogram                    
#hist(list,bins=5,density=True,color='b',label='data')
#show()
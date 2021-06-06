#importing libraries

import numpy as np
# for arrays
import matplotlib.pyplot as plt
# plotting
import pandas as pd
#dealing with data sets

######################################

#importing data sets

data_set1= pd.read_csv('x.csv')
x=data_set1.iloc[:,:].values
data_set2= pd.read_csv('y.csv')
y=data_set2.iloc[:,:].values
#print(x)
#print(y)

zeros=[]
x_filtered=[]
y_filtered=[]
for i in range(0,len(x)):
    n=list(x[i][2:]).count(0)
    m=list(y[i][2:]).count(0)
    if (n<=25 and m<=25):
        x_filtered.append(x[i][:])
        y_filtered.append(y[i][:])
    else:
        zeros.append(x[i][0])
#####################################

def myFunc(e):
  return e['r']
def myFun(e):
  return e['name']
  # function for sorting
####################################
from scipy import stats
# import scipy for coll.
c=[]
d=[]
# define dic. c for savim=ng names with CC
from numpy import nan
# define nan value
import math
# import math to detect nan values
from scipy.stats import spearmanr
const_var=[]
const_var2=[]
for i in range(0, len(x)):
    m=stats.pearsonr(x_filtered[i,2:],y_filtered[i,2:])
    n= spearmanr(x_filtered[i,2:],y_filtered[i,2:])
    if math.isnan(m[0]):
      const_var.append(x_filtered[i,0])
    else:
      c.append({'name':x_filtered[i,0],'r':m[0]})
    if math.isnan(n[0]):
      const_var2.append(x[i,0])
    else:
      d.append({'name':x_filtered[i,0],'r':n[0]})
# for loop produces CC for Dataset
c.sort(key=myFunc)
d.sort(key=myFunc)
print(c[0])
print(c[-1])
print(d[0])
print(d[-1])
# print max and min CC

for i in range(0, len(x)):
  if x_filtered[i,0]==myFun(c[0]):
    a=i
  if x_filtered[i,0]==myFun(c[-1]):
    b=i
for i in range(0, len(x)):
  if x_filtered[i,0]==myFun(d[0]):
    e=i
  if x_filtered[i,0]==myFun(d[-1]):
    f=i
    # Get the number for the genes containing max and min CC
plt.scatter(x_filtered[a,2:],y_filtered[a,2:])
plt.scatter(x_filtered[b,2:],y_filtered[b,2:])
#plt.legend()
plt.show()
plt.scatter(x_filtered[e,2:],y_filtered[e,2:])
plt.show()
plt.scatter(x_filtered[f,2:],y_filtered[f,2:])
#plt.legend()
plt.show()
"""
plt.scatter(x[a,2:],y[a,2:])
#plt.plot(x[a,2:],y[a,2:])
plt.xlabel('normal '+ myFun(c[0]))
plt.ylabel('cancerous '+ myFun(c[0]))
plt.title('min -ve CC exp. lvl')
plt.show()
#plotting the exp.lvl
plt.scatter(x[b,2:],y[b,2:])
#plt.plot(x[b,2:],y[b,2:])
plt.xlabel('normal '+ myFun(c[-1]))
plt.ylabel('cancerous '+ myFun(c[-1]))
plt.title('max +ve CC exp. lvl')
plt.show()
print(a,b)
"""
"""
ind_p=[]
rel_p=[]
count=0
import scipy

for i in range(0, len(x)):
  w=scipy.stats.ttest_ind(x[i,2:],y[i,2:],equal_var=True,nan_policy='propagate')
  ind_p.append(w[1])
#print(ind_p)

pr=[]
namer=[]
for i in range(0, len(x)):
  w=scipy.stats.ttest_rel(x[i,2:],y[i,2:],nan_policy='omit')
  if math.isnan(w[1]):
    print('y')
  else:
    #rel_p.append({'name':x[i,0],'r':w[1]})
    pr.append(w[1])
    namer.append(x[i,0])
    count=count+1
#for i in range(0, count):
#  if math.isnan(myFunc(rel_p[i])):
 #    # const_var.append(x[i,0])
  #   print('n')
    #else:
     # c.append({'name':x[i,0],'r':m[0]})
print(count)
#print(len(rel_p))
#print(myFunc(rel_p[]))
#print(pr)
import statsmodels.api
gene_names = ["G1", "G2", "G3", "G4", "G5"]
p = [0.1, 0.2, 0.01, 0.015, 0.04]
count_corr2=0
for i in range(0, len(p)):
  if p[i]<=0.05:
    #print(namer[i])
    count_corr2=count_corr2+1
print(count_corr2)
count_corr1=0
for i in range(0, len(pr)):
  if pr[i]<=0.05:
    #print(namer[i])
    count_corr1=count_corr1+1
print(count_corr1)
cr=statsmodels.stats.multitest.multipletests(pr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
#ci=statsmodels.stats.multitest.multipletests(ind_p, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
#print(cr[1][5])
count_corr=0
for i in range(0, len(pr)):
  if cr[1][i]<=0.05:
    #print(namer[i])
    count_corr=count_corr+1
print(count_corr)
#for i in range(0, len(x)):
 # if cr[i]<=0.025:
  #  print(x[i,0])

p = [0.1, 0.2, 0.01, 0.015, 0.04]
q_fdr =statsmodels.stats.multitest.multipletests(p, method = 'fdr_bh')
qq=q_fdr[1]
print(qq[2])
"""
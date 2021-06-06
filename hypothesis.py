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

from scipy import stats
# import scipy for coll.
from numpy import nan
# define nan value
import math
# import math to detect nan values

import scipy
pr=[]
namer=[]
no_pv_gener=[]
for i in range(0, len(x)):
  w=scipy.stats.ttest_rel(x[i,2:],y[i,2:],nan_policy='omit')
  if math.isnan(w[1]):
    no_pv_gener.append(x[i,0])
  else:
    pr.append(w[1])
    namer.append(x[i,0])
import statsmodels.api
gene_names = ["G1", "G2", "G3", "G4", "G5"]
p = [0.1, 0.2, 0.01, 0.015, 0.04]
count_corr1=0
name_deg1=[]
for i in range(0, len(pr)):
  if pr[i]<=0.05:
    name_deg1.append(namer[i])
    count_corr1=count_corr1+1
cr=statsmodels.stats.multitest.multipletests(pr, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
count_corr2=0
name_deg2=[]
for i in range(0, len(pr)):
  if cr[1][i]<=0.05:
    name_deg2.append(namer[i])
    count_corr2=count_corr2+1
count_tf=0
count_ft=0
deg_dif_ind=[]
for i in range(0,count_corr1):
    if name_deg1[i] not in name_deg2:
        count_tf=count_tf+1
        deg_dif_ind.append(name_deg1[i])
for i in range(0,count_corr2):
    if name_deg2[i] not in name_deg1:
        count_ft=count_ft+1
print(count_tf)
print(count_ft)
print(len(name_deg1))
print(len(name_deg2))
print(len(deg_dif_ind))
print(deg_dif_ind)
###########################################

import statsmodels.api
pi=[]
namei=[]
no_pv_genei=[]
for i in range(0, len(x)):
  w=scipy.stats.ttest_ind(x[i,2:],y[i,2:],nan_policy='omit')
  if math.isnan(w[1]):
    no_pv_genei.append(x[i,0])
  else:
    pi.append(w[1])
    namei.append(x[i,0])
gene_names = ["G1", "G2", "G3", "G4", "G5"]
p = [0.1, 0.2, 0.01, 0.015, 0.04]
count_corr3=0
name_deg3=[]
for i in range(0, len(pi)):
  if pi[i]<=0.05:
    name_deg3.append(namei[i])
    count_corr3=count_corr3+1
ci=statsmodels.stats.multitest.multipletests(pi, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
count_corr4=0
name_deg4=[]
for i in range(0, len(pi)):
  if ci[1][i]<=0.05:
    name_deg4.append(namei[i])
    count_corr4=count_corr4+1
count_tf2=0
count_ft2=0
deg_dif_rel=[]
for i in range(0,count_corr3):
    if name_deg3[i] not in name_deg4:
        deg_dif_rel.append(name_deg3[i])
        count_tf2=count_tf2+1
for i in range(0,count_corr4):
    if name_deg4[i] not in name_deg3:
        count_ft2=count_ft2+1
print(count_tf2)
print(count_ft2)
print(len(name_deg3))
print(len(name_deg4))
print(len(deg_dif_rel))
print(deg_dif_rel)

#######################################################
#point 4
name_common=[]
name_dist_rel=[]
name_dist_ind=[]
count_common=0
count_dist_rel=0
count_dist_ind=0
for i in range(0,len(name_deg2)):
    if name_deg2[i] in name_deg4:
        name_common.append(name_deg2[i])
        count_common=count_common+1
for i in range(0,len(name_deg2)):
    if name_deg2[i] not in name_common:
        name_dist_rel.append(name_deg2[i])
        count_dist_rel=count_dist_rel+1
for i in range(0,len(name_deg4)):
    if name_deg4[i] not in name_common:
        name_dist_ind.append(name_deg4[i])
        count_dist_ind=count_dist_ind+1
print(count_dist_ind)
print(count_dist_rel)
print(count_common)
print(name_dist_ind)
print(name_dist_rel)
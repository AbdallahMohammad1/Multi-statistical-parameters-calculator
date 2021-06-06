#importing libraries

import numpy as np
# for arrays
import matplotlib.pyplot as plt
# plotting
import pandas as pd
#dealing with data sets
import statsmodels.stats.multitest as multi

######################################

#importing data sets

data_set1= pd.read_csv('x.csv')
x=data_set1.iloc[:,:].values
data_set2= pd.read_csv('y.csv')
y=data_set2.iloc[:,:].values
#print(x)
#print(y)

zeros=[]
x_filt=[]
y_filt=[]
for i in range(0,len(x)):
    n=list(x[i][2:]).count(0)
    m=list(y[i][2:]).count(0)
    if (n<=25 and m<=25):
        x_filt.append(list(x[i][:]))
        y_filt.append(list(y[i][:]))
    else:
        zeros.append(x[i][0])
#print(x_filt[0][:])
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
# define dic. c for saving names with CC
from numpy import nan
# define nan value
import math
# import math to detect nan values
const_var=[]
for i in range(0, len(x_filt)):
    m=stats.pearsonr(x_filt[i][2:],y_filt[i][2:])
    if math.isnan(m[0]):
      const_var.append(x_filt[i][0])
    else:
      c.append({'name':x_filt[i][0],'r':m[0]})
# for loop produces CC for Dataset
c.sort(key=myFunc)
print(c[0])
print(c[-1])
# print max and min CC

for i in range(0, len(x_filt)):
  if x_filt[i][0]==myFun(c[0]):
    a=i
  if x_filt[i][0]==myFun(c[-1]):
    b=i

    # Get the number for the genes containing max and min CC
plt.scatter(x_filt[a][2:],y_filt[a][2:])
#plt.plot(x[a,2:],y[a,2:])
plt.xlabel('normal '+ myFun(c[0]))
plt.ylabel('cancerous '+ myFun(c[0]))
plt.title('min -ve CC exp. lvl')
plt.show()
#plotting the exp.lvl
plt.scatter(x_filt[b][2:],y_filt[b][2:])
#plt.plot(x[b,2:],y[b,2:])
plt.xlabel('normal '+ myFun(c[-1]))
plt.ylabel('cancerous '+ myFun(c[-1]))
plt.title('max +ve CC exp. lvl')
plt.show()
#print(a,b)
################################################### Independent ##########################################
p_val_ind_no_correction = []
nan_ind_no_correction = []
name_ind =[]
for i in range(0, len(x_filt)):
  f = stats.ttest_ind(x_filt[i][2:],y_filt[i][2:],nan_policy='propagate')
  if math.isnan(f[1]):
    nan_ind_no_correction.append(x_filt[i][0])
  else:
    p_val_ind_no_correction.append(f[1])
    name_ind.append(x_filt[i][0])

sig_ind_no_correction = 0
name_sig_ind_no_correction = []
for i in range(0, len(p_val_ind_no_correction)):
  if  p_val_ind_no_correction[i] <= 0.05:
    sig_ind_no_correction = sig_ind_no_correction + 1
    name_sig_ind_no_correction.append(name_ind[i])

q = multi.multipletests(p_val_ind_no_correction, method='fdr_bh')
p_val_ind_corrected = q[1]

sig_ind_corrected = 0
name_sig_ind_corrected = []
for i in range(0, len(p_val_ind_corrected)):
  if  p_val_ind_corrected[i] <= 0.05:
    sig_ind_corrected = sig_ind_corrected + 1
    name_sig_ind_corrected.append(name_ind[i])

number_ind_difference = sig_ind_no_correction - sig_ind_corrected
name_ind_difference = []
for i in range (0, sig_ind_no_correction):
  if name_sig_ind_no_correction[i] not in name_sig_ind_corrected:
    name_ind_difference.append(name_sig_ind_no_correction[i])


print('p_val_ind_number = ',len(p_val_ind_no_correction))
print('sig_ind_no_correction_number = ',sig_ind_no_correction)
print('sig_ind_corrected_number = ',sig_ind_corrected)
print('number_ind_difference = ',number_ind_difference)
#print('name_ind_difference = ',name_ind_difference)

###################################### Related ##############################################
p_val_rel_no_correction = []
nan_rel_no_correction = []
name_rel =[]
for i in range(0, len(x_filt)):
  f = stats.ttest_rel(x_filt[i][2:],y_filt[i][2:],nan_policy='propagate')
  if math.isnan(f[1]):
    nan_rel_no_correction.append(x_filt[i][0])
  else:
    p_val_rel_no_correction.append(f[1])
    name_rel.append(x_filt[i][0])

sig_rel_no_correction = 0
name_sig_rel_no_correction = []
for i in range(0, len(p_val_rel_no_correction)):
  if  p_val_rel_no_correction[i] <= 0.05:
    sig_rel_no_correction = sig_rel_no_correction + 1
    name_sig_rel_no_correction.append(name_rel[i])

q = multi.multipletests(p_val_rel_no_correction, method='fdr_bh')
p_val_rel_corrected = q[1]

sig_rel_corrected = 0
name_sig_rel_corrected = []
for i in range(0, len(p_val_rel_corrected)):
  if  p_val_rel_corrected[i] <= 0.05:
    sig_rel_corrected = sig_rel_corrected + 1
    name_sig_rel_corrected.append(name_rel[i])

number_rel_difference = sig_rel_no_correction - sig_rel_corrected
name_rel_difference = []
for i in range (0, sig_rel_no_correction):
  if name_sig_rel_no_correction[i] not in name_sig_rel_corrected:
    name_rel_difference.append(name_sig_rel_no_correction[i])

print('p_val_rel_number = ',len(p_val_rel_no_correction))
print('sig_rel_no_correction_number = ',sig_rel_no_correction)
print('sig_rel_corrected_number = ',sig_rel_corrected)
print('number_rel_difference = ',number_rel_difference)
#print('name_rel_difference = ',name_rel_difference)

####################### 4th requirement in hypothesis ########################################

number_sig_common_corrected = 0
name_sig_common_corrected = []
for i in range (0, sig_ind_corrected):
  if name_sig_ind_corrected[i] in name_sig_rel_corrected:
    number_sig_common_corrected = number_sig_common_corrected + 1
    name_sig_common_corrected.append(name_sig_ind_corrected[i])

number_sig_distinct_ind_corrected = 0
name_sig_distinct_ind_corrected = []
for i in range (0, sig_ind_corrected):
  if name_sig_ind_corrected[i] not in name_sig_common_corrected:
    number_sig_distinct_ind_corrected = number_sig_distinct_ind_corrected + 1
    name_sig_distinct_ind_corrected.append(name_sig_ind_corrected[i])

number_sig_distinct_rel_corrected = 0
name_sig_distinct_rel_corrected = []
for i in range (0, sig_rel_corrected):
  if name_sig_rel_corrected[i] not in name_sig_common_corrected:
    number_sig_distinct_rel_corrected = number_sig_distinct_rel_corrected + 1
    name_sig_distinct_rel_corrected.append(name_sig_rel_corrected[i])

print('number_sig_common_corrected = ',number_sig_common_corrected)
# print('name_sig_common_corrected = ',name_sig_common_corrected)
print('number_sig_distinct_ind_corrected = ',number_sig_distinct_ind_corrected)
#print('name_sig_distinct_ind_corrected = ',name_sig_distinct_ind_corrected)
print('number_sig_distinct_rel_corrected = ',number_sig_distinct_rel_corrected)
#print('name_sig_distinct_rel_corrected = ',name_sig_distinct_rel_corrected)


# p = [0.1,0.2,0.01,0.015,0.04]
# test = multi.multipletests(p,method='fdr_bh')
# print(test[1])

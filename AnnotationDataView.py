# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 13:29:29 2019

@author: CLARG38
"""

import sys
filein="C:/Users/CLARG38/OneDrive - The Toronto-Dominion Bank/Documents/Python Scripts/ParsedData_2018.xlsx"
from statsmodels.stats.outliers_influence import variance_inflation_factor
import os
import re
from shutil import copyfile
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler,OneHotEncoder, LabelEncoder
from datetime import datetime
import time
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import scipy.stats as st


prefix="C:/Users/CLARG38/Downloads/PlayDate/"
infile=pd.read_excel(os.path.join(prefix,"Processed","AnnotationInfoOnly_pared.xlsx"))
print infile.columns
infile=infile[(infile['Micro-injection aborted'] > 0) | (infile['Genotype confirmed'] > 0)]
infile=infile[~((infile['Micro-injection aborted'] > 0) & (infile['Genotype confirmed'] > 0))]

#infile=infile[(infile['Micro-injection aborted'] ==0 )]
#print infile.columns
#infile=infile[(infile['Micro-injection aborted'] > 0) & (infile['Genotype confirmed'] > 0)]

def backtoone(x):
    if x >= 1:
        return 1
    else:
        return 0
import scipy.stats as st
from sklearn.metrics import confusion_matrix    
    
import matplotlib.pyplot as plt

#for itm in ['minEntropy','degree_centrality','degree','GCcontent','Num_Window','Length','PercentageCpG','max_size','CpGsites','total_bp','Num_lowEntropy_Window','min_size','PolII','H3K4me1','num_annotated']:
#
#    infile=pd.read_excel(os.path.join(prefix,"Processed","AnnotationInfoOnly_pared.xlsx"))
#    
#    infile=infile[(infile['Micro-injection aborted'] > 0) | (infile['Genotype confirmed'] > 0)]
#    infile=infile[~((infile['Micro-injection aborted'] > 0) & (infile['Genotype confirmed'] > 0))]
#    #print infile
#    infile['Genotype confirmed']=infile['Genotype confirmed'].apply(lambda x: backtoone(x))
#    infile.drop('Micro-injection aborted',axis=1,inplace=True)
#    mylist=[itm,'Genotype confirmed']
#    infile=infile[mylist].copy()
#
#    infile.dropna(how='any',inplace=True)
#    gtC=infile[(infile['Genotype confirmed'] == 1)]
#    gtA=infile[(infile['Genotype confirmed'] == 0)]#len(infile)-gtC
#    if not len(infile):
#        continue
#    nA=round(np.mean(gtC[itm].values),3)
#    nB=round(np.mean(gtA[itm].values),3)
#    print len(gtC),len(gtA)
#    sA=round(np.std(gtC[itm].values),3)
#    #nA=round(np.mean(filter(lambda j: not np.isnan(j),gtC[itm].values),3))
#    #nB=round(np.mean(filter(lambda j: not np.isnan(j),gtA[itm].values),3))
#    
#    #sA=round(np.std(filter(lambda j: not np.isnan(j),gtC[itm].values),3))
#    
#    
#    target=infile[['Genotype confirmed']].copy().values.ravel()
# 
#    infile.drop('Genotype confirmed',axis=1,inplace=True)
#      
#    
#
#    model = RandomForestClassifier(n_estimators=1500, n_jobs=-1, random_state=412,max_depth=1)
#    model.fit(infile, target)
#    predicted = model.predict(infile)
#    
#    conf_mat = confusion_matrix(target, predicted)
#    print itm,len(infile),conf_mat,"\t",conf_mat[0:,0][0],len(gtA),"\t",conf_mat[0:,1][1],len(gtC)
#    #print gtC[itm].values
#    plt.hist(gtA[itm].values,bins=25)
#    plt.hist(gtC[itm].values,bins=15)
#    plt.show()
#    print nA,sA,nB
#    print "\n\n"
#sys.exit()

infile=pd.read_excel(os.path.join(prefix,"Processed","AnnotationInfoOnly_pared.xlsx"))

infile=infile[(infile['Micro-injection aborted'] > 0) | (infile['Genotype confirmed'] > 0)]
infile=infile[~((infile['Micro-injection aborted'] > 0) & (infile['Genotype confirmed'] > 0))]

infile['Genotype confirmed']=infile['Genotype confirmed'].apply(lambda x: backtoone(x))
infile.drop('Micro-injection aborted',axis=1,inplace=True)
infile.dropna(how='any',inplace=True)
rawtar=infile[['Genotype confirmed']].copy().values
target=infile[['Genotype confirmed']].copy().values.ravel()
 
infile.drop('Genotype confirmed',axis=1,inplace=True)

    


model = RandomForestClassifier(n_estimators=1500, n_jobs=-1, random_state=412,min_impurity_decrease=0.001)
model.fit(infile, target)
predicted = model.predict(infile)

conf_mat = confusion_matrix(target, predicted)
print conf_mat
print conf_mat[0:,0][0],conf_mat[0:,1][1]



importances = model.feature_importances_
indices = np.argsort(importances)
#indices.reverse()
#arr[::-1]
for j in indices[::-1]:
    print infile.columns[j], "=", importances[j]
##print list(targetvals)
#sys.exit()
sys.exit()
features_names = merged.columns
svm = svm.SVC(kernel='linear')
print svm,"...now fitting it"
svm.fit(merged, target)
f_importances(svm.coef_, features_names)
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

from ModelDistributions import FitDist
fitter=FitDist()


        
if os.path.exists(filein):
    print "we got it"
  
openexcel = pd.ExcelFile(filein)
dfs=openexcel.parse()   ##with empty argument just reads the first sheet?

dfs=dfs[dfs['CutSize'] < 5000]

cols=dfs.columns

def countGRNA(col):
    col=str(col)
    num=col.count(",")
    return num

dfs['numGRNA']=pd.Series(dfs['CutPositions'].map(lambda x: countGRNA(x)),index=dfs.index)

def binaryconvert(col):
    col=str(col)
    regtrue=re.compile("t|TRUE$|y|yes",flags=re.I|re.X)
    if regtrue.search(col):
        return 1
    else:
        return 0

def dateconvert(col):
    ##mostly format "yyyy-mm-dd"
    if len(str(col).strip()):
        dateNM=datetime.strptime(col, '%Y-%m-%d')
        timestamp=time.mktime(dateNM.timetuple())

    else:
        timestamp=0
    return timestamp


###Lets predict specific NaN columns using SVM classifier
 #                 'Allele Type',
  #                'Allele Subtype',

# =============================================================================
categorical=dfs[['Consortium_x',
                  'Production Centre_x',
                  'Allele Type',
                  'Allele Subtype',
                  'Status Name',
                  'Assay Carried Out',
                  'mRNA Nuclease',
                  'Protein Nuclease',
                  'Viability',
                  'numGRNA',
                  'Delivery Method'
                  ]].copy()

categorical=dfs[['Status Name',
                  'mRNA Nuclease',
                  'Protein Nuclease',
                  'Viability',
                  'numGRNA',
                  'Delivery Method'
                  ]].copy()
# =============================================================================
binary=dfs[['Experimental',
                 'Truncated_guide'
                 ]].copy()
countdata=dfs[['Voltage',
                'Protein Concentration',
                'mRNA Concentration',
                '#Embryos Injected',
                '#Embryos Survived',
                '#Embryos Transfered',
                'CutSize'
                 ]].copy()

datedata=dfs[['Mi Date'
        ]].copy()


###
count_col=countdata.columns
cat_col=categorical.columns
bin_col=binary.columns
date_col=datedata.columns


def explore_distribution(countdata):
    import matplotlib.pyplot as plt
    for c in count_col:
    
        condition=map(lambda x: x > 0 and x < 5000,countdata[c])
        ndat=np.extract(condition,countdata[c].values)      ##select only data that matches
                                                            ##above condition
        ##scale the number of bins to the data
        numbins=min(len(set(ndat)),max(len(set(ndat)),len(set(ndat))/10))
    
        best_fit_name, best_fit_params = fitter.best_fit_distribution(ndat, numbins)
        best_dist = getattr(st, best_fit_name)
        pdf=fitter.make_pdf(best_dist,best_fit_params)
        print c,"\t",best_fit_name, best_dist,best_fit_params


countdata.fillna(0,inplace=True)
categorical.fillna('Unknown',inplace=True)
binary=binary.applymap(lambda x: binaryconvert(x))
datedata=datedata.applymap(lambda x: dateconvert(x))


hdat=[]
for c in count_col:
    ndat=countdata[c]
    numbins=max(4,min(len(set(ndat)),min(len(set(ndat))/30,10)))
    bns=np.linspace(1,max(ndat),numbins)
    bns=np.insert(bns,0,0)  ##first bin is a catch-all for non-annotated/NaN values
    hdat.append(np.digitize(countdata[c], bns))

newcat=pd.DataFrame(np.vstack(hdat).T,columns=count_col,index=countdata.index)
#print newcat
#sys.exit()
def encodeCategorical(categorical,encoder="OneHot"):
    ##Pre-processing of categorical data
    cat_col=categorical.columns
    if encoder == "nominal":
        trans=[]
        for j in cat_col:
            le=LabelEncoder()
        
            le.fit(categorical[j].unique())
            transformed=le.transform(categorical[j].values)
            trans.append(transformed)
            #sys.exit()
        categorical=pd.DataFrame(zip(*trans),columns=cat_col,index=categorical.index)
    elif encoder == "OneHot":
        for j in cat_col:
            categorical = pd.concat([categorical,pd.get_dummies(categorical[j], prefix=j)],axis=1)
            categorical.drop([j],axis=1, inplace=True)
    return categorical

newcat=encodeCategorical(newcat)
categorical=encodeCategorical(categorical)
#print newcat.columns
#sys.exit()


#bins = np.linspace(1, 5000, 10)
#digitized = np.digitize(countdata["CutSize"], bins)
#print digitized
#for q,a in zip(countdata["CutSize"],digitized):
#    print q,a




scaled_counts = StandardScaler().fit_transform(countdata.values)
counts = pd.DataFrame(scaled_counts, index=countdata.index, columns=countdata.columns)



##Re-merge into a single data frame
#merged=pd.concat([categorical,binary,countdata,datedata],axis=1)

##in this case we have a binary matrix
merged=pd.concat([categorical,binary,newcat],axis=1)

from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics import mutual_info_score

mi_score = lambda a, b: mutual_info_score(a, b)


mergedC=merged.columns
for j in range(merged.shape[1]):
    for k in range(j+1,merged.shape[1]):
        datJ,datK=merged.values[j,:],merged.values[k,:]
        lclScore=adjusted_mutual_info_score(datJ,datK).round(3)
        cor,pval=st.spearmanr(datJ,datK)
        if lclScore > 0.2 and cor < 0:
            print mergedC[j],mergedC[k],lclScore,st.spearmanr(datJ,datK)
  

sys.exit()
vif = pd.DataFrame()
vif["VIF Factor"] = [variance_inflation_factor(countdata.values, i) for i in range(countdata.shape[1])]
vif["features"] = countdata.columns
print vif.round(1)
sys.exit()

def f_importances(coef, names):
    imp = coef
    imp,names = zip(*sorted(zip(imp,names)))
    for j,k in zip(imp,names):
        print k,":",j


##Prepare for training
#data=merged.copy()
#target=merged['Delivery Method_Unknown'].copy()
target=merged['Status Name_Genotype confirmed'].copy()

#queryFeature=['Delivery Method_Electroporation',
#              'Delivery Method_Cytoplasmic Injection',
#              'Delivery Method_Pronuclear Injection']
queryFeature=['Status Name_Genotype confirmed','Status Name_Micro-injection aborted','Status Name_Founder obtained']
for q in queryFeature:
    merged.drop([q],axis=1,inplace=True)
##data.drop(['Delivery Method_Unknown'],axis=1,inplace=True)
#categorical.drop(['Delivery Method_Unknown'],axis=1,inplace=True)

#print target

#targetvals=target.values#.ravel()
bayesTarget=dfs[['Delivery Method']].copy()
bayesTarget.fillna('Unknown',inplace=True)
def MNB(data,target):
    #print data
    from sklearn.naive_bayes import BernoulliNB,MultinomialNB,GaussianNB,ComplementNB
    from sklearn import metrics
    #gnb = MultinomialNB()
    #gnb = ComplementNB()
    gnb = BernoulliNB(alpha=0.0,fit_prior=False)
    #gnb = GaussianNB()
    model = gnb.fit(data.values, target)
    print dir(model)

    #exponent=np.e
    exponent=np.e
    #model_prior=map(lambda l: np.exp(l),model.class_log_prior_)
    #model_prior=map(lambda l: np.power(exponent,l),model.class_log_prior_)
    model_prior=model.class_log_prior_
    #print model.class_log_prior_
    #print model_prior
    
    predicted = model.predict(data.values)
    
    pos_prob=model.feature_log_prob_[0, :]
    neg_prob=model.feature_log_prob_[1, :]
    

    matrix=data.values
    #print list(matrix[0,:]) == list(matrix[0])
    #print list(matrix[1,:]) == list(matrix[1])
    #print matrix[:,0]
    #print matrix[:,1]


    pos_prob_l=map(lambda l: np.power(exponent,l), model.feature_log_prob_[0, :])
    neg_prob_l=map(lambda l: np.power(exponent,l), model.feature_log_prob_[1, :]) 
    

    #sumprob=
    indices_pos = np.argsort(pos_prob)
    indices_neg = np.argsort(neg_prob) 

    totalUS=0
    totalPR=0
    for idx,pred in zip(range(len(data)),predicted):
        posvec=merged.iloc[idx].values
        #negvec=merged.iloc[idx].map(lambda x: invert(x)).values

        sum_pos=pos_prob*posvec#/float(model_prior[1])#*float(model_prior[1])
        sum_neg=neg_prob*posvec#/float(model_prior[0])#*float(model_prior[0])
        #if sum_pos+sum_neg < 8.52:
        #    print pred
        #    totalUS+=1
        
        ##pseudocount
        prodP=sum(sum_pos)+model_prior[1]
        prodN=sum(sum_neg)+model_prior[0]

        if prodN > prodP:
            totalUS+=1
            print "\t",idx,pred,"T",prodP,prodN,prodP > prodN,abs(prodN-prodP)

        totalPR+=pred
    print totalUS
    print totalPR

    #print metrics.classification_report(target, predicted)
    #print metrics.confusion_matrix(target, predicted)


def invert(x):
    return abs(1-x)
#MNB(merged,target)


clf = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial').fit(merged, target)
#print dir(clf)
for j,k in zip(merged.columns,clf.coef_[0]):
    print j,float(k)
sys.exit()


model = RandomForestClassifier(n_estimators=1500, n_jobs=-1, random_state=42)
model.fit(merged, target)
importances = model.feature_importances_
indices = np.argsort(importances)
#indices.reverse()
#arr[::-1]
for j in indices[::-1]:
    print merged.columns[j], "=", importances[j]
##print list(targetvals)
#sys.exit()
sys.exit()
features_names = merged.columns
svm = svm.SVC(kernel='linear')
print svm,"...now fitting it"
svm.fit(merged, target)
f_importances(svm.coef_, features_names)
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



# =============================================================================
#             Mi Attempt URL
#             Mi Attempt External Ref
#             Gene Marker Symbol_x
#             Gene MGI Accession ID_x
#             Consortium_x
#             Production Centre_x
#             Status Name
#             Report Micro Injection Progress To Public
#             IS Active
#             Experimental
#             Zygote
#             IMPC Mutant Strain Donor
#             MF External Ref
#             Plasmid Generated Guides
#             gRNA Concentrations
#             gRNA Concentrations Individually Set?
#             gRNA Sequence (+ strand)
#             Chromosome (+ strand)
#             Start Co-od
#             End Co-od
#             Truncated_guide
#             Individually set gRNA concentrations
#             gRNA Sequence (+ strand).1
#             Chromosome (+ strand).1
#             Start Co-od.1
#             End Co-od.1
#             Truncated_guide.1
#             Individually set gRNA concentrations.1
#             gRNA Sequence (+ strand).2
#             Chromosome (+ strand).2
#             Start Co-od.2
#             End Co-od.2
#             Truncated_guide.2
#             Individually set gRNA concentrations.2
#             gRNA Sequence (+ strand).3
#             Chromosome (+ strand).3
#             Start Co-od.3
#             End Co-od.3
#             Truncated_guide.3
#             Individually set gRNA concentrations.3
#             mRNA Nuclease
#             mRNA Concentration
#             Protein Nuclease
#             Protein Concentration
#             Delivery Method
#             Voltage
#             #Pulses
#             Mi Date
#             #Embryos Injected
#             #Embryos Survived
#             Embryo Transfer Day
#             #Embryos Survived to 2 cell stage
#             #Embryos Transfered
#             #Founder Pups Born
#             #Founders Assayed
#             Assay Carried Out
#             #Founders Selected For Breeding
#             #G0 with detected mutation
#             #G0 NHEJ event detected
#             #G0 deletion event detected
#             #G0 HR event detected
#             #G0 HDR event detected
#             #G0 all donor insertions detected
#             #G0 subset of donors inserted detected
#             Mi Attempt URL.1
#             Gene Marker Symbol_y
#             Gene MGI Accession ID_y
#             Consortium_y
#             Production Centre_y
#             F1 Colony Name
#             Genotype Confirmed
#             Released from Genotyping (WTSI)
#             Report F1 Colony To Public
#             Background Strain
#             Uploaded Trace file
#             MGI Allele Symbol Superscript
#             MGI Allele Accession ID
#             MGI Identifier
#             Allele Type
#             Allele Subtype
#             Sequences
#             CutPositions
#             CutSize
#             #G0 with INFERRED mutation
#             #Embryos Transfered INFERRED
#             #Pups Born INFERRED
#             Viability
#             FounderRate_transf
#             GLTRate
#             numGRNA
# =============================================================================





if os.path.exists(filein):
    print "we got it"
  
openexcel = pd.ExcelFile(filein)
dfs=openexcel.parse()   ##with empty argument just reads the first sheet?

def countGRNA(col):
    col=str(col)
    num=col.count(",")
    return num

def binaryconvert(col):
    col=str(col)
    regtrue=re.compile("t|TRUE$|y|yes",flags=re.I|re.X)
    if regtrue.search(col):
        return 1
    else:
        return 0




def unique_E(arr):
    return len(filter(lambda x: x == "Genotype confirmed",arr)),len(arr)
def first(arr):
    return list(arr)[0]
    
dfs['numGRNA']=pd.Series(dfs['CutPositions'].map(lambda x: countGRNA(x)),index=dfs.index)

if __name__ == "__main__":
    

    cols=dfs.columns

    grouped = dfs.groupby('Gene Marker Symbol_x')
    groupeddf=grouped.agg({
            'Gene MGI Accession ID_x':first,
            '#Embryos Injected': sum,
            '#Embryos Survived': sum,
            '#Embryos Survived to 2 cell stage': sum,
            '#Embryos Transfered': sum,
            '#Founder Pups Born': sum,
            '#Founders Assayed': sum,
            '#Founders Selected For Breeding': sum,
            '#G0 with detected mutation': sum,
            '#G0 NHEJ event detected': sum,
            '#G0 deletion event detected': sum,
            '#G0 HR event detected': sum,
            '#G0 HDR event detected': sum,
            '#G0 all donor insertions detected': sum,
            '#G0 subset of donors inserted detected': sum,
            'Status Name':unique_E,
            'Viability':first
            })
    print groupeddf['Viability']
    vib=groupeddf[groupeddf['Viability'].isin(['Viable','Subviable','Lethal'])]
    subset=vib[['Gene MGI Accession ID_x','Viability']].copy()
    subset.to_excel("C:/Users/CLARG38/Downloads/PlayDate/Viability.xlsx")
#            'number A': 'sum',
#            'number B': 'min'})






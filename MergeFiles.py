# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:05:33 2019

@author: CLARG38
"""
import string, re,sys
import os
from collections import defaultdict
import glob
from cPickle import dump,load
import pandas as pd
import platform
import numpy as np
if platform.node().startswith("greg-clarks-macbook"):
    prefix="/Users/clarkgr1/Desktop/aws"

        
if __name__ == "__main__":
    
    ##If we add a list of genes, we compile for only these
    #
    main=pd.read_excel(os.path.join(prefix,"GeneInformation.xlsx"))
    symbols=list(main['Symbol'].values)
    newannot=[]
    with open(os.path.join(prefix,"check_MGI.csv")) as io:
        for line in io:
            data=map(lambda s: s.strip("\""),line.strip().split(","))
            uid,s,na=data
            if s in symbols:
                lclDF=main[main['Symbol']==s].copy()
                if lclDF['DEG_essentiality'].values[0] in [0,1]:
                    newannot.append(lclDF.values[0][0])
    subset=main[main['MGI_ID'].isin(newannot)].copy()
    subset.to_excel('Subset.xlsx',index=False)


        

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


from ParseFiles import Annotation,Gene



def grab_sequences(excelfile):
    filedf=pd.read_excel(excelfile)
    ##get list of genes - convert to string from numpy returned unicode
    symbols=list(map(lambda s: str(s),filedf['Gene Marker Symbol_x'].values))
    return symbols


from itertools import islice
def window(seq, n=200):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def percentGC(seq):
    gc=seq.count('G')+seq.count('C')
    try:
        content=gc/float(len(seq))
    except ZeroDivisionError:
        return 0
    return content


def CpGisland(seq):
    lclFC=percentGC(seq)
    observed=seq.count('CG')
    c=seq.count('C')
    g=seq.count('G')
    try:
        expected=((c+g)/2.)**2/200.
        ratio=observed/float(expected)
    except ZeroDivisionError:
        ratio=0  
    if lclFC > 0.5 and ratio > 0.6:
        return 1
    else:
        return 0
import math
def entropy(string):
        # get probability of chars in string
        prob = [ float(string.count(c)) / len(string) for c in dict.fromkeys(list(string)) ]
        # calculate the entropy
        entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob ])
        return entropy    
#from scipy.stats import entropy

if __name__ == "__main__":
    
    ##If we add a list of genes, we compile for only these
    #
    inputfile="C:/Users/CLARG38/Downloads/PlayDate/Aggregated.xlsx"
    genelist=grab_sequences(inputfile)
    #print genelist[:20]
    lcl=Annotation()

    
    symbol2mgi=load(open(os.path.join(lcl.prefix,"Processed",'symbol2mgi.pkl'),'rb'))
    uid2mgi=load(open(os.path.join(lcl.prefix,"Processed",'uid2mgi.pkl'),'rb'))
    protein2mgi=load(open(os.path.join(lcl.prefix,"Processed",'protein2mgi.pkl'),'rb'))
    ensembl2mgi=load(open(os.path.join(lcl.prefix,"Processed",'ensembl2mgi.pkl'),'rb'))
    
#    genes=load(open(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_GeneInit.pkl"),'rb'))
#    regulatory=load(open(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_data_regulatory.pkl"),'rb'))
#    print "read regulatory"
#    degree=load(open(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_geneDegree.pkl"),'rb'))
#
    sequences=load(open(os.path.join(lcl.prefix,"Processed","Sequences.pkl"),'rb'))
#
#    lncRNA=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"lncrna_.xlsx"))
#    enhancers=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_Enhancers.xlsx"))
#    promoters=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_Promoters.xlsx"))
#    miRNA=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_mirna.xlsx"))
    
    genedata=pd.read_excel(os.path.join(lcl.prefix,"Processed","_Sequences.xlsx"))
    
    import numpy as np
    inputfile="C:/Users/CLARG38/Downloads/PlayDate/Aggregated.xlsx"
    maindf=pd.read_excel(inputfile)
    seqinfo={}
    capture=[]
    cols=['MGI_ID','Length','GCcontent','CpGsites','PercentageCpG','Num_Window','minEntropy','Num_lowEntropy_Window']
    cnt=0
    for symbol,seqdata in sequences.iteritems():
        try:
            mgi=symbol2mgi[symbol]
        except KeyError:
            try:
                mgi=ensembl2mgi[seqdata[2]]
            except KeyError:
                pass
        try: 
            mgi
        except NameError: 
            continue
        sequence=seqdata[1]
        if not len(sequence):
            continue
        splitLarge=["".join(x) for x in window(sequence, 200)]
        CpGislands=map(lambda l: CpGisland(l),splitLarge)
        Einfo=map(lambda l: entropy(l),splitLarge)
        try:
            CpGcoverage=sum(CpGislands)/float(len(CpGislands))
        except ZeroDivisionError:
            CpGcoverage=0
        try:
            minE=min(Einfo)
        except ValueError:
            minE=1.92
        LowEntropy=filter(lambda e: e < 1,Einfo)
        GCcontent=percentGC(sequence)
        print cnt,"\t",mgi,round(CpGcoverage,3),np.std(Einfo),minE,len(LowEntropy)
        capture.append([mgi,len(sequence),GCcontent,sum(CpGislands),CpGcoverage,len(CpGislands),minE,len(LowEntropy)])    
        cnt+=1
    #print capture
    print len(capture[0])
    df=pd.DataFrame(capture,columns=cols)
    df.to_excel(os.path.join(lcl.prefix,"Processed","_Sequences.xlsx"))    

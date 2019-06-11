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
    symbols=list(map(lambda s: str(s),filedf['Gene MGI Accession ID_x'].values))
    return symbols


if __name__ == "__main__":
    
    ##If we add a list of genes, we compile for only these
    #
    inputfile="C:/Users/CLARG38/Downloads/PlayDate/Aggregated.xlsx"
    genelist=grab_sequences(inputfile)
    #print genelist[:20]
    lcl=Annotation()

    def convert(symbol):
        try:
            mgi=symbol2mgi[symbol]
        except KeyError:
            mgi=None
        return mgi

 
    ##Supporting data
    symbol2mgi=load(open(os.path.join(lcl.prefix,"Processed",'symbol2mgi.pkl'),'rb'))
    uid2mgi=load(open(os.path.join(lcl.prefix,"Processed",'uid2mgi.pkl'),'rb'))
    protein2mgi=load(open(os.path.join(lcl.prefix,"Processed",'protein2mgi.pkl'),'rb'))
    ensembl2mgi=load(open(os.path.join(lcl.prefix,"Processed",'ensembl2mgi.pkl'),'rb'))
    
    genes=load(open(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_GeneInit.pkl"),'rb'))
    
    ####
    for g,d in genes.iteritems():
        print g,dir(d)
    sys.exit()
    
    regulatory=pd.read_excel(os.path.join(lcl.prefix,"Processed","Regulatory.xlsx"))

    #df = df.merge(df_, on='join_col_name')

    degreeINFO=load(open(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_geneDegree.pkl"),'rb'))
    degdata=[]
    degcolumns=['MGI_ID','degree_centrality','degree']
    for i,k in degreeINFO.iteritems():
        degdata.append([i]+k)
    degree=pd.DataFrame(degdata,columns=degcolumns)
    print len(regulatory)
    master=pd.merge(regulatory,degree,how='outer',on=['MGI_ID'])
    print len(degree)
    print len(master)
    lncRNA=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"lncrna_.xlsx"))

    lncRNA.loc[:,'brain'] = lncRNA.loc[:,'brain_1'].add(lncRNA.loc[:,'brain_2'])
    lncRNA.drop(['brain_1','brain_2'],axis=1,inplace=True)
 
    lncRNA.loc[:,'heart'] = lncRNA.loc[:,'heart_1'].add(lncRNA.loc[:,'heart_2'])
    lncRNA.drop(['heart_1','heart_2'],axis=1,inplace=True)
 
    lncRNA.loc[:,'liver'] = lncRNA.loc[:,'liver_1'].add(lncRNA.loc[:,'liver_2'])
    lncRNA.drop(['liver_1','liver_2'],axis=1,inplace=True)

#    u'MGI_ID', u'chr', u'start', u'end', u'name', u'strand', u'exonNo',
#       u'exonSize', u'exonStart', u'brain_rep1_exp_fpkm',
#       u'brain_rep2_exp_fpkm', u'heart_rep1_exp_fpkm', u'heart_rep_exp_fpkm2',
#       u'intestine_rep1_exp_fpkm', u'intestine_rep_exp_fpkm2',
#       u'kidney_rep1_exp_fpkm', u'kidney_rep2_exp_fpkm',
#       u'liver_rep1_exp_fpkm', u'liver_rep2_exp_fpkm', u'spleen_rep1_exp_fpkm',
#       u'spleen_rep2_exp_fpkm', u'testes_rep1_exp_fpkm',
#       u'testes_rep2_exp_fpkm', u'thymus_rep1_exp_fpkm',
#       u'thymus_rep2_exp_fpkm', u'es_rep1_exp_fpkm', u'es_rep2_exp_fpkm',
#       u'esg_rep1_exp_fpkm', u'esg_rep2_exp_fpkm', u'esg_rep3_exp_fpkm',
#       u'length', u'RNA_Size', u'ORF_Size', u'Ficket_Score', u'Hexamer_Score',
#       u'CPAT_Score', u'PhyloCSF_Score', u'annotation',
#       u'annotatio_intergenic

    lncRNA=lncRNA[['MGI_ID','brain','heart','liver','Ficket_Score','ORF_Size']].copy()
    Lagg=lncRNA.groupby('MGI_ID')
    Lseries=Lagg.agg({'MGI_ID':'first',
            'brain': 'mean',
            'heart': 'mean',
            'liver': 'mean',
            'Ficket_Score': 'mean',
            'ORF_Size': 'mean'
            })
    Lcols=Lseries.columns
    lncRNA=pd.DataFrame(Lseries,columns=Lcols)
    

    master=pd.merge(master,lncRNA,how='outer',on=['MGI_ID'])
    print len(master)
    enhancers=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_Enhancers.xlsx"))
    enhancers['MGI_ID']=enhancers['gene'].apply(lambda x: convert(x))       ##convert symbol to mgi id
    enhancers = enhancers[pd.notnull(enhancers['MGI_ID'])]      ##Remove rows where no annotation
    enhancers['E14.5-brain'].fillna(0,inplace=True)
    enhancers['E14.5-heart'].fillna(0,inplace=True)
    enhancers['E14.5-limb'].fillna(0,inplace=True)
    Eagg=enhancers.groupby('MGI_ID')
    
    ##u'chromosome', u'enhancer', u'TSS', u'gene', u'transcript',
    ##   u'E14.5-brain', u'E14.5-heart', u'E14.5-limb
    Eseries=Eagg.agg({
            'MGI_ID':'first',
            'E14.5-brain': max,
            'E14.5-heart': max,
            'E14.5-limb': max,
            'chromosome': 'first',
            'enhancer': 'first',
            'TSS': 'first',
            'transcript': 'first'
            })
    cols=Eseries.columns
    enhancers=pd.DataFrame(Eseries,columns=cols)
    enhancers=enhancers[['MGI_ID','E14.5-brain','E14.5-heart','E14.5-limb']].copy()
    
    master=pd.merge(master,enhancers,how='outer',on=['MGI_ID'])
    
    
    promoters=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_Promoters.xlsx"))
    promoters['MGI_ID']=promoters['Gene'].apply(lambda x: convert(x))       ##convert symbol to mgi id
    promoters = promoters[pd.notnull(promoters['MGI_ID'])]      ##Remove rows where no annotation
    #u'Gene', u'Transcript',        u'chr',  u'direction',
     #         u'TSS',    u'E-brain',    u'E-heart',     u'E-limb',
     #      u'MGI_ID
    promoters['E-brain'].fillna(0,inplace=True)
    promoters['E-heart'].fillna(0,inplace=True)
    promoters['E-limb'].fillna(0,inplace=True)
    Pagg=promoters.groupby('MGI_ID')
    Pseries=Pagg.agg({
        'MGI_ID':'first',
        'E-brain': max,
        'E-heart': max,
        'E-limb': max,
        'chr': 'first',
        'direction': 'first',
        'TSS': 'first',
        'Transcript': 'first'
        })
    Pcols=Pseries.columns
    promoters=pd.DataFrame(Pseries,columns=Pcols)
    promoters=promoters[['MGI_ID','E-brain','E-heart','E-limb']].copy()
    
    master=pd.merge(master,promoters,how='outer',on=['MGI_ID'])
    
    miRNA=pd.read_excel(os.path.join(lcl.prefix,"Processed",lcl.fileprefix+"_mirna.xlsx"))
    miRNA.rename(columns={'mgiID':'MGI_ID'}, inplace=True)

    master=pd.merge(master,miRNA,how='outer',on=['MGI_ID'])
    
    seqdata=pd.read_excel(os.path.join(lcl.prefix,"Processed","_Sequences.xlsx"))
    

    master=pd.merge(master,seqdata,how='outer',on=['MGI_ID'])
    inputfile="C:/Users/CLARG38/Downloads/PlayDate/Aggregated.xlsx"
    maindf=pd.read_excel(inputfile)
    maindf.rename(columns={'Gene MGI Accession ID_x':'MGI_ID'}, inplace=True)
    print maindf.columns
    print len(maindf),len(master)
    maindf=pd.merge(maindf,master,how='left',on=['MGI_ID'])
    #maindf.dropna(subset=['MGI_ID'],axis='rows')
    maindf.to_excel(os.path.join(lcl.prefix,"Processed","AnnotationInfo.xlsx"),index=False)
    print len(maindf)
    print maindf
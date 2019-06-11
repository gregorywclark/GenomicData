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


class Annotation():
    """The class servers as a container for some properties and will be 
    the parent node and container for all Gene instances"""


    if platform.node().startswith("greg-clarks-macbook"):
        prefix="/Users/clarkgr1/Desktop/aws"
    elif platform.node().startswith("WIN"):
        prefix="C:/Users/CLARG38/Downloads/PlayDate/"

    def __init__(self,genelist=[]):
        self.genelist=genelist
        self.definedgenes=bool(len(genelist))
        self.fileprefix='_partial' if len(genelist) > 0 else '_all'
        self.chromosomes=list(map(lambda s: str(s),range(1,20))) + ['X','Y','MT']
        self.genes={}
        self.uid2mgi={}
        self.symbol2mgi={}
        self.protein2mgi={}
        self.ensembl2mgi={}

        self.reannotate=False

        self.link_ids_=False
        self.miRNA_=False
        self.lcRNA_=False
        self.encode_=False
        self.viability_=False       #functional and viability info
        self.ppi_=False

        self.regflag=False
        self.seqflag=False
        self.OOGE_added=False
        ##Collections from excel spreadsheets are put into data-frames
        ##Will be reconciled/normalized later
        self.lncRNA=pd.DataFrame()

    def link_ids(self):
        filename=os.path.join(self.prefix,'Processed',self.fileprefix+'_GeneInit.pkl')
        if os.path.exists(filename) and not self.reannotate:
            print "Gene info exists, retrieving it..."
            self.genes=load(open(filename,'rb'))
            if len(self.genes) < len(self.genelist):
                print "We need to add more"
            else:
                self.link_ids_=True
                return
            
        ##We have self.genes to act as a container for Gene instances

        markers=os.path.join(self.prefix,"MRK_Sequence.rpt")
        ncbiInfo=os.path.join(self.prefix,"MGI_Gene_Model_Coord.rpt")

        refprot=re.compile("ENSMUSP[0-9]{3,15}",flags=re.I|re.X)
        ensgene=re.compile("ENSMUSG[0-9]{3,15}",flags=re.I|re.X)
        with open(markers) as mrk:
            mrk.next()
            #skip the header
            for line in mrk:
                data=line.split("\t")
                mgi,symbol,refseq=data[0],data[1],data[-4]
                ##Master Annotation
                self.symbol2mgi[symbol]=mgi
                
                refseq=filter(lambda l: refprot.search(l),refseq.split("|"))
                if self.definedgenes:
                    if mgi in self.genelist or symbol in self.genelist:
                        ins=Gene(symbol)
                        ins.mgi=mgi
                        ins.symbol=symbol
                        ins.proteins=refseq
                        for protein in refseq:
                            self.protein2mgi[protein]=mgi
                        self.genes[mgi]=ins
                else:
                    ##This means doing ALL genes 
                    ins=Gene(symbol)
                    ins.mgi=mgi
                    ins.symbol=symbol
                    ins.proteins=refseq
                    for protein in refseq:
                        self.protein2mgi[protein]=mgi
                    self.genes[mgi]=ins

        with open(ncbiInfo) as ncbi:
            ncbi.next()
            ##skip the header
            for line in ncbi:
                data=line.strip().split("\t")
                mgi,symbol,uid,ens,chrom,start,end,strand=map(lambda j: data[j],[0,2,5,10,11,12,13,14])

                ##We must have defined mgi via markers, or skip iteration
                try:    self.genes[mgi]
                except KeyError:    continue
                
                ##Master Annotation
                self.uid2mgi[uid]=mgi
                if ensgene.search(ens):
                    self.ensembl2mgi[ens]=mgi
                    self.genes[mgi].ensembl=ens
                    
                self.genes[mgi].uid=uid
                self.genes[mgi].chromosome=chrom
                

                try:
                    self.genes[mgi].chromosome=chrom
                    self.genes[mgi].start=int(start)
                    self.genes[mgi].end=int(end)
                    self.genes[mgi].length=abs(int(end)-int(start))+1
                    self.genes[mgi].strand=strand

                except ValueError:
                    try:
                        self.genes[mgi].chromosome=data[6]
                        self.genes[mgi].start=int(data[7])
                        self.genes[mgi].end=int(data[8])
                        self.genes[mgi].length=abs(int(data[8])-int(data[7]))+1
                        self.genes[mgi].strand=data[9]
                        self.genes[mgi].co_ord='ensembl'

                    except ValueError:
                        ##Strike abhorrent abberation
                        self.genes[mgi].co_ord=None
 
        self.genes={i:j for i,j in self.genes.iteritems() if j.chromosome and j.co_ord}
        
        if not self.definedgenes:
            ioA=open(os.path.join(self.prefix,"Processed",'symbol2mgi.pkl'),'wb')
            ioB=open(os.path.join(self.prefix,"Processed",'uid2mgi.pkl'),'wb')
            ioC=open(os.path.join(self.prefix,"Processed",'protein2mgi.pkl'),'wb')
            ioD=open(os.path.join(self.prefix,"Processed",'ensembl2mgi.pkl'),'wb')
            dump(self.symbol2mgi,ioA)
            dump(self.uid2mgi,ioB)
            dump(self.protein2mgi,ioC)
            dump(self.ensembl2mgi,ioD)
            ioA.close()
            ioB.close()
            ioC.close()
            ioD.close()
            
        filehandler=open(filename,'wb')
        dump(self.genes,filehandler)
        filehandler.close()
        self.link_ids_=True


    def add_viability(self):

        assert self.link_ids_, "Must link ids prior to annotation at gene level"
        ##Functional annotation from ncbi
        with open(os.path.join(self.prefix,"generifs_basic")) as rif:
            for line in rif:
                data=line.strip().split("\t")
                if data[0] == "10090":  ## ncbi taxonomy id for mouse 
                    if re.search("lethal|fatal|dead|viable|essential",data[-1]) and not re.search("protect|resistance|rescue",data[-1]):
                        try:
                            ##We have to convert from uid to mgi (our chosen ID) first
                            mtch=re.search("lethal|fatal|dead|viable|essential",data[-1])
                            #print mtch.group(0),data[-1].split()
                            ant=mtch.group(0)
                            if not ant:
                                print "\n\n\n\n"
                                print re.search("lethal|fatal|dead|viable|essential",data[-1])
                            lclmgi=self.uid2mgi[data[1].strip()]
                            self.genes[lclmgi].annotation=ant
                        except KeyError:
                            pass

        ##Viability information from impc               
        with open(os.path.join(self.prefix,"Viability.csv")) as vb:
            for line in vb:
                data=line.strip().split(",")
                if data[-1] not in ['Viable','Subviable','Lethal']:
                    continue
                try: 
                    self.genes[data[-2]].viability=data[-1]
                except KeyError: 
                    pass
        self.viability_=True

    def add_OOGE(self):
        self.OOGE_added=True
        with open(os.path.join(self.prefix,"Mus_musculus_OOGE.csv")) as deg:
            header=deg.next().strip().split(",")
            ##locus,symbols,datasets,datasetIDs,essentiality status,essentiality consensus
            for line in deg:
                data=line.strip().split(",")
                ensb,symbol,gclass=data[0],data[1],data[-2]
                       
                try:
                    mgi=self.symbol2mgi[symbol]
                    
                except KeyError:
                    try:
                        mgi=self.ensembl2mgi[ensb]
                    except KeyError:
                        mgi=None
                    
                if not mgi:
                    continue
                try:
                    self.genes[mgi]
                except KeyError:
                    continue
                if gclass == "NE":       
                    self.genes[mgi].OOGEclass=0
                else:
                    self.genes[mgi].OOGEclass=1
                    
                
    def add_DEG(self):
        assert self.OOGE_added,"Must run MM to augment gene identifier relationships"
        ##
        with open(os.path.join(self.prefix,"DEG_essential_mouse.csv")) as deg:
            ##no header here
            for line in deg:
                data=line.strip().split("\t")
                ensb,gclass=data[3:5]
                try:
                    mgi=self.ensembl2mgi[ensb]
                    self.genes[mgi]
                    if gclass == "NE":       
                        self.genes[mgi].DEGclass=0
                    else:
                        self.genes[mgi].DEGclass=1
                except KeyError:
                    continue
    def add_OOGE_human(self):
        assert self.OOGE_added,"Must run MM to augment gene identifier relationships"
        ##locus,symbols,datasets,datasetIDs,essentiality status,essentiality consensus
        def clean_class(split):
            line=split.strip("\"")
            cleanline=line.split(",")
            classes=list(set(cleanline))
            if len(set(classes)) == 1:
                if classes[0] == "NE":
                    return 0
                elif classes[0] == "E":
                    return 1
            else:
                return classes.count("E")/float(len(classes))
            
        with open(os.path.join(self.prefix,"Homo_sapiens_OOGE.csv")) as deg:
            ##no header here
            for line in deg:
                data=line.strip().split(",")
                #print data
                symbol,gclass=data[1],data[-2]
                symbol=symbol.capitalize()
                try:
                    mgi=self.symbol2mgi[symbol]
                    if data[-1].strip() == "Nonessential":
                        gclass=0
                    elif data[-1].strip() == "Essential":
                        gclass=1
                    else:
                        dinfo=line.split("\"")[-2].split(",")
                        gclass=round(dinfo.count("E")/float(len(dinfo)),3)
                    self.genes[mgi].HS_OOGEclass=gclass
                    #print symbol,mgi,self.genes[mgi].HS_OOGEclass
                except KeyError:
                    continue
        
        
class Gene(Annotation):
    
    """Creating a data construct that will capture all information about a gene.
    However we take a protein-centric view because there are multiple proteins found 
    for each gene identifier"""

    refprot=re.compile("ENSMUSP[0-9]{3,15}",flags=re.I|re.X)

    def __init__(self,symbol):
        ##ids obtain from ncbi and mgi
        self.id=symbol
        self.proteins=[]
        self.mgi=None
        self.symbol=None
        self.uid=None
        self.ensembl=None
        
        ##ncbi co-ordinates
        self.co_ord='ncbi'
        self.chromosome=0
        self.start=0
        self.end=0
        self.strand=''
        
        ##metrics based on protein network
        self.degree=None
        self.dc=None
        
        self.annotation=None
        self.viability=None
        

        self.regulatory=[]
        self.DEGclass=None
        self.OOGEclass=None
        self.HS_OOGEclass=None

def grab_sequences(excelfile):
    filedf=pd.read_excel(excelfile)
    ##get list of genes - convert to string from numpy returned unicode
    symbols=list(map(lambda s: str(s),filedf['Gene Marker Symbol_x'].values))
    return symbols
        
if __name__ == "__main__":
    
    ##If we add a list of genes, we compile for only these
    #
    inputfile="C:/Users/CLARG38/Downloads/PlayDate/Aggregated.xlsx"
    genelist=grab_sequences(inputfile)


    ann=Annotation(genelist)
    ann.fileprefix='essential'
    ann.link_ids()
    ann.add_OOGE()
    ann.add_viability()
    ann.add_DEG()
    ann.add_OOGE_human()
    
    vb=[]
    lt=[]
    data=[]
    col=['MGI_ID','Symbol','EnsemblID','Viability','DEG_essentiality','Human_Essentiality(OOGE)','GeneRif_keyword']
    for gene,i in ann.genes.iteritems():
        if not i.viability and (i.OOGEclass or i.HS_OOGEclass):  
            print gene,i.symbol,i.ensembl,i.viability,i.annotation,i.DEGclass,i.OOGEclass,i.HS_OOGEclass
        data.append([gene,i.symbol,i.ensembl,i.viability,i.DEGclass,i.HS_OOGEclass,i.annotation])
    df=pd.DataFrame(data,columns=col)
    df.to_excel(os.path.join(ann.prefix,'Processed','GeneInformation.xlsx'),index=False)

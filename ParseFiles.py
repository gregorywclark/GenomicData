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
        self.fileprefix='_partial' if len(genelist) > 0 else ''
        self.chromosomes=list(map(lambda s: str(s),range(1,20))) + ['X','Y','MT']
        self.genes={}
        self.uid2mgi={}
        self.symbol2mgi={}

        self.reannotate=False

        self.link_ids_=False
        self.miRNA_=False
        self.lcRNA_=False
        self.encode_=False
        self.viability_=False       #functional and viability info
        self.ppi_=False

        self.regflag=False
        self.seqflag=False
        ##Collections from excel spreadsheets are put into data-frames
        ##Will be reconciled/normalized later
        self.lncRNA=pd.DataFrame()

    def link_ids(self):
        filename=os.path.join(self.prefix,'Processed',self.fileprefix+'_GeneInit.pkl')
        if os.path.exists(filename) and not self.reannotate:
            print "Gene info exists, retrieving it..."
            self.genes=load(open(filename,'rb'))
            self.link_ids_=True
            print "done."

        ##We have self.genes to act as a container for Gene instances

        markers=os.path.join(self.prefix,"MRK_Sequence.rpt")
        ncbiInfo=os.path.join(self.prefix,"MGI_Gene_Model_Coord.rpt")

        refprot=re.compile("ENSMUSP[0-9]{3,15}",flags=re.I|re.X)
        refgene=re.compile("ENSMUSG[0-9]{3,15}",flags=re.I|re.X)
        with open(markers) as mrk:
            mrk.next()
            #skip the header
            for line in mrk:
                data=line.split("\t")
                mgi,symbol,refseq=data[0],data[1],data[-4]
                self.symbol2mgi[symbol]=mgi
                refseq=filter(lambda l: refprot.search(l),refseq.split("|"))
                if self.definedgenes:
                    if mgi in self.genelist:
                        ins=Gene(symbol)
                        ins.mgi=mgi
                        ins.symbol=symbol
                        ins.proteins=refseq
                        self.genes[mgi]=ins
                else:

                        ins=Gene(symbol)
                        ins.mgi=mgi
                        ins.symbol=symbol
                        ins.proteins=refseq
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
                

                self.uid2mgi[uid]=mgi
                self.genes[mgi].uid=uid
                self.genes[mgi].chromosome=chrom
                self.genes[mgi].ensembl=ens
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
                    if re.search("lethal|fatal|dead|viable",data[-1])and not re.search("protect|resistance|rescue",data[-1]):
                        try:
                            ##We have to convert from uid to mgi (our chosen ID) first
                            lclmgi=self.uid2mgi[data[1].strip()]
                            self.genes[lclmgi].annotation=data[-1]
                        except KeyError:
                            pass

        ##Viability information from impc               
        with open(os.path.join(self.prefix,"Viability.csv")) as vb:
            for line in vb:
                data=line.strip().split(",")
                if data[-1] not in ['Viable','Subviable','Lethal']:
                    continue
                try: 
                    self.genes[data[-2]]=data[-1]
                except KeyError: 
                    pass
        self.viability_=True

    def add_lncrna(self):

        assert self.link_ids_, "Must link ids prior to annotation at gene level"
        import pandas as pd
        
        filein=os.path.join(self.prefix,"Processed",self.fileprefix+"lncrna_.xlsx")
        if os.path.exists(filein):
            print "We've already processed lncRNA...retriving"
            self.lncRNA=pd.read_excel(filein)
            print "Got it"
            return
        
        ###File related to:
        ##Chromatin and RNA Maps Reveal Regulatory Long Noncoding RNAs in Mouse
        ##Gireesh K. Bogu, Pedro Vizán, Lawrence W. Stanton, Miguel Beato, Luciano Di Croce, Marc A. Marti-Renom
        ##Mol Cell Biol. 2015 Dec 28;36(5):809-19. 
        ##doi: 10.1128/MCB.00955-15.

        lncfile=os.path.join(self.prefix,"zmb999101146so2_M.xls")      
        #chrmtn_lncfile=os.path.join(self.prefix,"zmb999101146so2.xls")

        ###columns for lncfile .... 
        ##chr	start	end   name	strand	exonNo	exonSize	exonStart
        ##brain_rep1_exp_fpkm	brain_rep2_exp_fpkm	heart_rep1_exp_fpkm	heart_rep_exp_fpkm2	intestine_rep1_exp_fpkm	intestine_rep_exp_fpkm2	
        ##kidney_rep1_exp_fpkm	kidney_rep2_exp_fpkm	liver_rep1_exp_fpkm	liver_rep2_exp_fpkm	spleen_rep1_exp_fpkm	spleen_rep2_exp_fpkm	
        ##testes_rep1_exp_fpkm	testes_rep2_exp_fpkm	thymus_rep1_exp_fpkm	thymus_rep2_exp_fpkm	es_rep1_exp_fpkm	es_rep2_exp_fpkm	
        ##esg_rep1_exp_fpkm	esg_rep2_exp_fpkm	esg_rep3_exp_fpkm	name	
        ##length	RNA_Size  ORF_Size	Ficket_Score	Hexamer_Score	CPAT_Score	Coding_Label	PhyloCSF_Score	annotation	 annotatio_intergenic

        datacollection=[]
        lncdf=pd.ExcelFile(lncfile)
        for sheet in lncdf.sheet_names:
            sheetdf=pd.read_excel(lncfile, sheet_name=sheet)
            sheet_cols=['MGI_ID']+list(sheetdf.columns)
            for idx,row in sheetdf.iterrows():
                ##column names in excel sheet for location are 'chr','start','end'
                ##Here we are iterating over lnRNA, then filtering our genes to find if there is
                ##an overlap between the two. If overlap, place data into a dataframe associated Annotation class
                chromosome,start,end=row['chr'][3:],int(row['start']),int(row['end'])
                lk={k: v for k, v in self.genes.iteritems() if v.chromosome == chromosome and v.start <= end and start < v.end}
                if len(lk):
                    for i,j in lk.iteritems():
                        datacollection.append([i]+list(row.values))
        self.lncRNA=pd.DataFrame(datacollection,columns=sheet_cols)
        self.lncRNA.to_excel(filein)
        
    def add_encode(self):
        assert self.link_ids_, "Must link ids prior to annotation at gene level"
    #### Files related to:
        #A map of the cis-regulatory sequences in the mouse genome
        #Yin Shen,1,* Feng Yue,1,* David F. McCleary,1 Zhen Ye,1 Lee Edsall,1 Samantha Kuan,1 Ulrich Wagner,1 Jesse Dixon,1,2,3 Leonard Lee,1 Victor V. Lobanenkov,4 and Bing Ren1,5
        #Nature. 2012 Aug 2; 488(7409): 116–120.
        #doi: 10.1038/nature11243
        
        ##the Bing Ren lab is a Major contributer to the ENCODE project. Majority of data for mm10 belongs to this lab
        ## We're using a processed version of this data from their paper
        ## An EPU is an enhancer-promoter unit     
        promoterUsage=os.path.join(self.prefix,"Ren_supplementary_table4_M.xlsx")        ##parse this as binary matrix gene promoter usage by each tissue type
        ##E-brain	  E-heart	E-limb

        promoterfile=os.path.join(self.prefix,"Processed",self.fileprefix+"_Promoters.xlsx")
        if os.path.exists(promoterfile):
            print "Promoters file exists...retrieving"
            prom_dfs=pd.read_excel(promoterfile)
            self.promoters=prom_dfs
            print "Got it"
        else:
            ##We could put these in a separate function
            ##But whatever
            def binary_func(x):
                if re.search("Yes|Y|True",x):   return 1
                else:   return 0
    
            promoters_in=pd.read_excel(promoterUsage)
            promoters=promoters_in[['Gene','Transcript','chr','direction','TSS','E-brain','E-heart','E-limb']].copy()
    
            apply_binary = ['E-brain','E-heart','E-limb']
            promoters[apply_binary] = promoters[apply_binary].applymap(lambda s: binary_func(s))
            promoters.to_excel(promoterfile,index=False)
            self.promoters=promoters
        
        enhancerfile=os.path.join(self.prefix,"Processed",self.fileprefix+"_Enhancers.xlsx")
        if os.path.exists(enhancerfile):
            print "Enhancers file exists...retrieving"
            sheet_dfs=pd.read_excel(enhancerfile)
            self.enhancers=sheet_dfs
            print "Got it"
            return
        enhancers=os.path.join(self.prefix,"Ren_supplementary_table7_M.xlsx")          ##parse this as binary matrix gene promoter usage by each tissue type
        embryo_sheets=["E14.5-brain","E14.5-heart","E14.5-limb"]
        for sheet in embryo_sheets:
            lclsheet=pd.read_excel(enhancers, sheet_name=sheet)
            lclsheet[sheet]=1       ##sheet is Tissue Name, assigned all genes in the sheet as found in this tissue
            try:
                ##merge all data into each-other, keep tissue
                sheet_dfs = pd.merge(sheet_dfs,lclsheet, how='outer', on=['enhancer', 'gene','TSS','chromosome','transcript'])          
            except NameError:
                sheet_dfs=lclsheet
  
        sheet_dfs.to_excel(enhancerfile,index=False)
        self.enhancers=sheet_dfs

    def add_miRNA(self):
        assert self.link_ids_, "Must link ids prior to annotation at gene level"
        ##From TargetScanMouse
        ##Most Mammalian mRNAs Are Conserved Targets of MicroRNAs 
        ##Robin C Friedman, Kyle Kai-How Farh, Christopher B Burge, David P Bartel.     
        ##Genome Research, 19:92-105 (2009). 
        site_scores=os.path.join(self.prefix,"Conserved_Site_Scores.txt")
    
    def add_regulatory(self):
        assert self.link_ids_, "Must link ids prior to annotation at gene level"
        ##from Ensembl BioMart
        from operator import itemgetter
        
        if self.regflag:
            return
        assert self.regflag == False,"We shouldn't be going this route otherwise"
        
        regulatory=os.path.join(self.prefix,"Regulatory_mart_export.csv")
        ##Chromosome        Start      End        Feature_type    Feature_type_class
      

        regulatorydata=os.path.join(self.prefix,"Processed",self.fileprefix+"_data_regulatory.pkl")
        container={c:[] for c in self.chromosomes}

        ###This is a large file ~7 million lines so we need to make some accomodations to parse
        ## segregate by chromosome,sort by start point, use interpolation to slice
        if os.path.exists(regulatorydata):
            print "Regulatory data file exists...retrieving it"
            regdata=open(regulatorydata,'rb')
            container=load(regdata)
        else:
            with open(regulatory) as rgl:
                header=rgl.readline().strip().split("\t")
                header.insert(0,"GeneID")
                for line in rgl:
                    data=line.strip().split()      
                    try:
                        chromosome,start,end=data[0].strip(),int(data[1]),int(data[2])
                        container[chromosome].append([start,end,data[3]])               
                    except (IndexError,KeyError):
                        ##Last line is Index Error (unknown why), multiple scaffolds are ignored (KeyError)
                        continue
            for c,v in container.iteritems():
                container[chromosome]=sorted(v,key=itemgetter(0))
                
            regdata=open(regulatorydata,'wb')
            dump(container,regdata)
            regdata.close()
                
        models={}      
        for c,v in container.iteritems():
            startrow=zip(*container[c])[0]
            def func(bp):
                return int((bp-min(startrow))*len(startrow)/float(max(startrow)-min(startrow)))
            models[c]=func

        completed=[]
        kk=0
        for g,gene in self.genes.iteritems():
            
            slcd=container[gene.chromosome]#[slicestart:sliceend]
            #print gene.mgi,"\t",gene.length,"\t",len(slcd),len(container[gene.chromosome]),"\t\t",slicestart,sliceend,"\t\t",gene.start,gene.end
            completed.append(gene.mgi)
            for s in slcd:
                if gene.start <= s[1] and s[0] < gene.end:
                    gene.regulatory.append(s)
            kk+=1
            if kk%1000 == 0:
                print kk,round(kk/len(self.genes),3)
        self.regflag=True        
        filename=os.path.join(self.prefix,"Processed",self.fileprefix+"_Genelist.pkl")
        filehandler=open(filename,'wb')
        dump(self.genes,filehandler)
        filehandler.close()
        
        
        #6829001 total rows


    def generate_sequences(self):
        assert self.link_ids_, "Must link ids prior to annotation at gene level"
        #assert self.genelist, "Must have a genelist defined, don't want this for all genes"
        assert self.seqflag == False, "Must not have completed already"
        
        import ast
        import requests
        
        seqfile=os.path.join(self.prefix,"Processed",self.fileprefix+"_Sequences.pkl")
        if os.path.exists(seqfile):
            self.sequences=load(open(seqfile),'rb')
            return
        server = "https://rest.ensembl.org"
        
        
        self.sequences={}
        cnt=1
        for k,g in self.genes.iteritems():
            gene=g
            symbol=k
            if gene.length > 50000 or str(gene.ensembl) == "null":
                continue
            print cnt,"\t",gene.mgi,gene.symbol, gene.length,gene.ensembl
            ext = "/lookup/id/%s?expand=1" % (str(gene.ensembl),)
            try:
                r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                if not r.ok:
                    r.raise_for_status()
            except:
                ##Will catch multiple errors
                print "ERROR"
                continue
    
            decoded = r.json()
            geneinfo=ast.literal_eval(repr(decoded))

            chromosome=geneinfo['seq_region_name']
            start=str(geneinfo['start'])
            end=str(geneinfo['end'])
            strand=str(geneinfo['strand'])
    
            SeqReq = "/sequence/region/mouse/%s:%s..%s:%s?content-type=text/x-fasta" %(chromosome,start,end,strand)
            #identifier=chromosome+":"+start+"-"+end+":"+strand
            try:
                r2 = requests.get(server+SeqReq, headers={ "Content-Type" : "text/x-fasta"})
                if not r2.ok:
                    r2.raise_for_status()
            except:
                continue
            sequence=str(r2.text).splitlines()
            sequence.pop(0)
            gene.sequence="".join(sequence)
            print gene.sequence[:25]
            self.sequences[gene.mgi]="".join(sequence)
            cnt+=1
        self.seqflag=True
        io=open(seqfile,'wb')
        dump(self.sequences,io)
        io.close()

        
    def add_ppi():

        MGInet=os.path.join(self.prefix,"MGI_ppi.tab")

        import networkx as nx
        G=nx.Graph()
        import numpy as np
        import scipy.stats as st
        degree=defaultdict(int)
        with open(MGInet) as mgi:
            for line in mgi:
                data=line.strip().split()
                G.add_edge(data[0],data[1])
                degree[data[0]]+=1
                degree[data[1]]+=1
        dc=nx.degree_centrality(G)
        Cv,Cl,Cs=[],[],[]
        for j,k in dc.iteritems():
            if j in viable:
                print j,k
                Cv.append(k)
            elif j in lethal:
                Cl.append(k)
            elif j in sub:
                Cs.append(k)

        ##As expected higher degree in Lethal > Subviable > Viable
        print np.mean(Cv),np.mean(Cs),np.mean(Cl)

        Cv.sort()
        Cl.sort()


class Gene(Annotation):
    
    """Creating a data construct that will capture all information about a gene.
    However we take a protein-centric view because there are multiple proteins found 
    for each gene identifier"""

    refprot=re.compile("ENSMUSP[0-9]{3,15}",flags=re.I|re.X)

    def __init__(self,symbol):
        ##ids obtain from ncbi and mgi
        self.id=symbol
        self.proteins=[]
        self.mgi=''
        self.symbol=''
        self.uid=''
        
        ##ncbi co-ordinates
        self.co_ord='ncbi'
        self.chromosome=0
        self.start=0
        self.end=0
        self.strand=''
        
        ##metrics based on protein network
        self.degree=0
        self.dc=0
        
        self.annotation=''
        self.viability=''
        

        self.regulatory=[]
        
        
        
if __name__ == "__main__":
    
    ##If we add a list of genes, we compile for only these
    #
    mylist=['MGI:1923631',
        'MGI:1921468',
        'MGI:1916649',
        'MGI:1925664',
        'MGI:1919828',
        'MGI:1914642',
        'MGI:1924311',
        'MGI:2149781',
        'MGI:1918645',
        'MGI:1924955']
    
    ann=Annotation()
    ann.link_ids()
#    ann.add_lncrna()
#    ann.add_encode()
    #ann.add_regulatory()
    ann.generate_sequences()

        

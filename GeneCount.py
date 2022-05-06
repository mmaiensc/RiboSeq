#!/usr/bin/env python
'''
Created on Jun 27, 2011
Finished on Jun 27, 2011
@author: eugeneoh
Genecount01 -> individualAnalyis
            -> normaliziationFactor
Determine number of reads that align per gene.
    [1] gene name; [2] strand; [3] low; [4] high; [5] raw counts; [6] rpM normalized \
    [7] rpkM normalized
    
2011.8.27 -> Add normFactor output. 
'''
import sys
class BaseRead:
    def __init__(self, base, read):
        self.base = base
        self.read = read
        
    def __str__(self):
        return '%d %d' % (self.base, self.read)
    
def uploadBaseRead(filenameForward,filenameReverse):
    '''Upload center-weighted base reads'''
    
    f=open(filenameForward,'r')
    forward={}
    for line in f:
        fields=line.split()
        #print "%s %s %s %d %f" %(fields[0], fields[1], fields[2], int(fields[0])+1, float(fields[1]))
        x=BaseRead(int(fields[0]), float(fields[1]))
        forward[x.base]=x.read
    f.close()
    
    r=open(filenameReverse,'r')
    reverse={}
    for line in r:
        fields=line.split()
        x=BaseRead(int(fields[0]), float(fields[1]))
        reverse[x.base]=x.read
    r.close()
    
    return forward,reverse
class GeneAnnotation:
    
    def __init__(self, name, low, high, strand, notes):
        self.name = name
        self.low = low
        self.high = high
        self.strand = strand
        self.notes = notes
        
    def __str__(self):
        return '%s %d %d %s %s' % (self.name, self.low, self.high, self.strand, self.notes)
def geneCount(filenameForward,filenameReverse,genelist):
    '''Count reads per gene here!
       Report total number of reads that align to coding sequences (per Million)'''
    
    (forward,reverse)=uploadBaseRead(filenameForward,filenameReverse)
    tscore=0
    l=[]
    f=open(genelist,'r')
    for line in f:
        fields=line.split('\t')
        if fields[3]=='-' and (fields[1]>fields[2]):
                temp_swap=fields[1]
                fields[1]=fields[2]
                fields[2]=temp_swap
        x=GeneAnnotation(fields[0], int(fields[1]), int(fields[2]), fields[3], \
                         fields[4].rstrip('\n'))
        if x.strand=='+':
            gscore=0
            for position in range(x.low,x.high+1):
                if position in forward: 
                    gscore+=forward[position]
                    tscore+=forward[position]
            l.append([x.name,x.strand,x.low,x.high,gscore,x.notes])
        elif x.strand=='-':
            gscore=0
            for position in range(x.low,x.high+1):
                if position in reverse:
                    gscore+=reverse[position]
                    tscore+=reverse[position]
            l.append([x.name,x.strand,x.low,x.high,gscore,x.notes])
    f.close()
    
    tscore=tscore/1e6
    return l,tscore
class GscoreVariables:
    
    def __init__(self, name, strand, low, high, gscore):
        self.name = name
        self.strand = strand        
        self.low = low
        self.high = high
        self.gscore = gscore
        
    def __str__(self):
        return '%s %s %d %d %d' % (self.name, self.strand, self.low, self.high, \
                                   self.gscore)        
        
def individualAnalysis(filenameForward,filenameReverse,genelist,output,normfactor):
    
    (l,tscore)=geneCount(filenameForward,filenameReverse,genelist)
    f=open(output,'w')
    for item in l:
        x=GscoreVariables(item[0],item[1],item[2],item[3],item[4])
        length=(x.high-x.low+1)/1e3
        f.write(x.name+'\t'+x.strand+'\t'+str(x.low)+'\t'+str(x.high)+'\t'+str(x.gscore)+'\t'+ \
                str(x.gscore/tscore)+'\t'+str(x.gscore/tscore/length)+'\n')
    f.close()
    
    f=open(normfactor,'w')
    f.write(str(tscore))
    f.close    
    
if __name__=='__main__':
    
    ### Change so it can take reference annotation and all output file names ###
    #genome='CDS_MG1655_mod_ann'    
    fnameForward=sys.argv[1]
    fnameReverse=sys.argv[2]
    genome_ann=sys.argv[3]
    output=sys.argv[4]
    output_norm=sys.argv[5]
    individualAnalysis(fnameForward,fnameReverse,genome_ann,output,output_norm)
   #  ### Create the following folders, if applicable ###
#     # Genecount01 -> individualAnalysis, normalizationFactor
# 
#     ### Do not change items below ###
#     path='/home/pkanabar/k/analysis/bowtie_mapping_MG1655/mapping/tertiary_analysis'
#     ### path='/home/pkanabar/k/analysis/bowtie_mapping_MG1655/mapping/trial_genecout01A'
#     glist='/home/pkanabar/k/reference/CDS_From_Eugene/'+genome
#     rawdata1=path+'/assignCount/'+readout
#         
#         
#         
#     ### Execute ###
#     output=path+'/Genecount01A/individualAnalysis/'+readout+'_IA'
#     norm=path+'/Genecount01A/normalizationFactor/'+readout+'_N'
#     individualAnalysis(fnameForward,fnameReverse,glist,output,norm)
    

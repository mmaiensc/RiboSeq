#!/usr/bin/env python
'''
Created on Apr 25, 2011
Finished on Apr 25, 2011
Edited on Apr 27, 2011
@author: eugeneoh
Count number of reads (weighted by the central residues) that have aligned to the genome.
Also count the number of these central residues each time for other analyses.
Determine 5' and 3' bias of sequence reads.
Determine length distribution of sequence reads.
Note: Rawdata00.py (2011.src) is computationally expensive [i.e. def GetAlignedReads()]!
      This revision substantially improves execution speed!
Modified Sep 20, 2018 by Mark Maienschein-Cline (mmaiensc@uic.edu), to allow user to give
minbound and maxbound params in command line.
'''
import sys
class AlignedReads:
    '''Define each aligned read to a class instance with various attributes'''
    
    def __init__(self, sequence, uniqueness, length, fiveprime, strand, mismatch):
        self.sequence = sequence
        self.uniqueness = uniqueness
        self.length = length
        self.fiveprime = fiveprime
        self.strand = strand
        self.mismatch = mismatch
    def __str__(self):
        return '%s %d %d %d %s %d' % (self.sequence, self.uniqueness, self.length, \
                                      self.fiveprime, self.strand, self.mismatch)
        
'''original counts: minbound=24,maxbound=46'''
def assignCount(distance,filename,minbound,maxbound):
    '''Count by 3prime minus distance'''
    f=open(filename,'r')
    forward={}
    reverse={}
    for line in f:
        fields=line.split()
        x=AlignedReads(fields[0], int(fields[1]), int(fields[2]), int(fields[3]), \
                       fields[4], int(fields[5]))
        if x.uniqueness==1 and minbound<x.length<maxbound and x.mismatch<3:
            if x.strand=='+':
                i=x.fiveprime+x.length-1-int(distance)
                if i in forward:
                	forward[i][0]+=1.0
                	forward[i][1]+=1
                else: 
                    forward[i]=[1.0,1] 
            elif x.strand=='-':
                j=x.fiveprime-x.length+1+int(distance)
                if j in reverse:
                	reverse[j][0]+=1.0
                	reverse[j][1]+=1
                else: 
                    reverse[j]=[1.0,1]
    f.close()
    
    return forward,reverse
def biasFrequency(filename,minbound,maxbound):
    '''Find nucleotide frequency at the 5' and 3' ends'''
    
    f=open(filename,'r')
    total = 0
    d5={'A':0,'C':0,'G':0,'T':0,'N':0}
    d3={'A':0,'C':0,'G':0,'T':0,'N':0}
    for line in f:
        fields=line.split()
        x=AlignedReads(fields[0], int(fields[1]), int(fields[2]), int(fields[3]), \
                       fields[4], int(fields[5]))
       # x=AlignedReads(fields[1], int(fields[3]), int(fields[4]), int(fields[7]), \
                      # fields[8], int(fields[10]))
        if x.uniqueness==1 and minbound<x.length<maxbound and x.mismatch<3:
            total+=1
            d5[x.sequence[0]]+=1
            d3[x.sequence[-1]]+=1
    #print sorted(d5.items()), sorted(d3.items())
    return total,d5,d3
            
def lengthDistribution(filename,minbound=0,maxbound=46):
    '''Determine length distribution'''
    
    f=open(filename,'r')
    d={}
    for line in f:
        fields=line.split()
        x=AlignedReads(fields[0], int(fields[1]), int(fields[2]), int(fields[3]), \
                       fields[4], int(fields[5]))
       # x=AlignedReads(fields[1], int(fields[3]), int(fields[4]), int(fields[7]), \
                #       fields[8], int(fields[10]))
        if x.uniqueness==1 and minbound<x.length<maxbound and x.mismatch<3:
            if x.length not in d: d[x.length]=1
            else: d[x.length]+=1
    #print sorted(d.items())
    return d    
    
def writeFiles(distance,filename,output1,output2,output3,output4,minbound,maxbound):
    '''Write physical text files
       output1=forward_strand; output2=reverse_strand; output3=bias_freq; output4=length_dist'''    
    (forward,reverse)=assignCount(distance,filename,minbound,maxbound)
    f=open(output1,'w')
    for l in sorted(forward.items()):
        f.write(str(l[0])+'\t'+str(l[1][0])+'\t'+str(l[1][1])+'\n')
    f.close()
    r=open(output2,'w')
    for l in sorted(reverse.items()):
        r.write(str(l[0])+'\t'+str(l[1][0])+'\t'+str(l[1][1])+'\n')
    r.close()
     
    (total,d5,d3)=biasFrequency(filename,minbound,maxbound)
    fn=open(output3,'w')
    for key in d5:
        fn.write(key+ \
                 "\t5': "+str('%0.3f'%(d5[key]/float(total)))+ \
                 "\t3': "+str('%0.3f'%(d3[key]/float(total)))+'\n')
    fn.write('Total reads: '+str(total))
    fn.close()
    
    d=lengthDistribution(filename)
    fn=open(output4,'w')
    for key in sorted(d.items()):
        fn.write(str(key[0])+'\t'+str(key[1])+'\n')
    fn.close()
    
if __name__ == '__main__':
	file=sys.argv[2]
	o1=sys.argv[3]
	o2=sys.argv[4]
	o3=sys.argv[5]
	o4=sys.argv[6]
	dist=sys.argv[1]
	minbound=int(sys.argv[7])
	maxbound=int(sys.argv[8])
	writeFiles(dist,file,o1,o2,o3,o4,minbound,maxbound)
	
    ### Change these variables as desired ###
    # date=sys.argv[1] 
#     path=sys.argv[2] 
#     file=sys.argv[3]
#     readout=sys.argv[4]
		
    ### Create the following folders ###
    # Rawdata01 -> assignCount, biasFrequency, lengthDistribution
    ### Do not change items below ### 
    #file='/Users/eugeneoh/Desktop/RawFiles/'+date+'/'+align
    #path='/Users/eugeneoh/Desktop/Solexa/'+date+'/Rawdata01/' 
    #o1=path+'assignCount/'+readout+'_F'
    #o2=path+'assignCount/'+readout+'_R'
    #o3=path+'biasFrequency/'+readout
    #o4=path+'lengthDistribution/'+readout
    #file='/home/pkanabar/projects/kannanK/try1/RawFiles/'+date+'/'+align
	
    
    #file='/home/pkanabar/projects/kannanK/Nov26_analysis/remove_known/tertiary_ana/trial_data'
    #file='/home/pkanabar/projects/kannanK/Nov26_analysis/remove_known/Matt_index6_LR_GT_ET_16_final_aligned.sam'
    #path='/home/pkanabar/projects/kannanK/Nov26_analysis/remove_known/tertiary_ana/'   
    # o1=path+'assignCount/'+readout+'_F'
#     o2=path+'assignCount/'+readout+'_R'
#     o3=path+'biasFrequency/'+readout
#     o4=path+'lengthDistribution/'+readout
 
    ### Execute ###
    
    

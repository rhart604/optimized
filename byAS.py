import sys,os
from Bio import SeqIO
import itertools
import datetime
import time
import argparse

importSAMPLE_ID = None
BARCODE_FILE = None
LOG_EVERY_N = 10000
SEL_STRING = 'paired' #global-perhaps make this an arg in future versions
VERSION = 'Version 3.0 - EDITED'
DESCRIPTION = '''\
  The goal is to separate single-cell RNAseq fastq files into two
  species specific subsets.  We iterate through the
  mixed-genome bam file to select read IDs matching the extracted barcodes
  along with alignment scores, choosing the best score by genome.
  Read IDs matching both genomes are deleted.
  Finally, open paired-end fastq files (must be in current directory)
  and output barcode-matching records
  to one genome, non-matching to an alternate genome. Requires that
  IDs match in paired files.
  '''

def checkOptions():

    parser = argparse.ArgumentParser(usage='samtools view <BAM> | splitByScore.py',description=DESCRIPTION)
    parser.add_argument('-v','--version', action='version', version=VERSION)
    parser.add_argument('-g','--genome',dest="sp1",default='hg38',
                    help='Selected genome identifier (default hg38)')
    parser.add_argument('-a','--altgenome',dest="sp2",default='mm10',
                    help='Non-matching genome identifier (default mm10)')
    parser.add_argument('-s','--sample',dest="sample_name",required=True,default=None,
                    help='Sample ID for insertion into file name')
    parser.add_argument('-e','--extension',dest='fastqext',default='fastq',
                    help='Fastq file extension (default fastq)')
    parser.add_argument('-f','--fastqdir',dest='fastqdir',default='fastq',
                    help='Directory with fastq files (default fastq)')
    parser.add_argument('-o','--output',dest="ofname",default=None,
                    help='File to log runtime output (default stdout)')

    args = parser.parse_args()
    return args

def getTag(tag, lis):
    '''
    Get optional tag from a SAM record.
    '''
    for t in lis:
        spl = t.split(':')
        if spl[0] == tag:
            return spl[-1]
    return None

def numFmt(value):
	return '{:,}'.format(value)

def fastqOut(genome,genDict,fastqDir,fastqRoot):
	fastqInName = fastqRoot 
	fastqOutName = fastqRoot.replace(SEL_STRING,genome) 
	fastq_iter = SeqIO.parse(open(fastqDir + '/' +fastqInName),'fastq')
	outFile = open(fastqDir + '/' + fastqOutName, 'w')
	ctr = 0
	foundCtr = 0
	startTime = time.time()
	for rec in fastq_iter:
			ctr += 1 
			if (rec.id in genDict): # record from chosen genotype
				foundCtr += 1
				SeqIO.write(rec,outFile,'fastq')
			else: # record not in dict
				continue
	fastq_iter.close()
	outFile.close()
	endTime = time.time() 
	print("Read " + numFmt(ctr) + " records from " + fastqInName + " and output " + numFmt(foundCtr) + " records to " + fastqOutName + " in " + '{:.2f}'.format(endTime - startTime) + " seconds.")
	return None

def main():
    '''Main.'''
    opt = checkOptions()
    sp1 = opt.sp1
    sp2 = opt.sp2
    sampleName = opt.sample_name
    fastqExt = opt.fastqext
    fastqDir = opt.fastqdir
    vString = 'v3.0 - HG/RH'
    oFname = opt.ofname
  
  #create lists
    readEnds = [ "_1", "_2" ] #perhaps add this as an option in the future
    genomes = [ sp1, sp2 ]

  # initialize counts
    ReadCount = 0 # total reads input
    GoodCount = 0 # total matching reads input
    UnmappedCount = 0
    MissingAS = 0
    MissingCB = 0 # reads without a CB tag
    sp1Count = 0 # reads assigned to sp1 
    sp2Count = 0 # reads assigned to sp2
  
  # create dictionaries
    sp1IDs = dict() #dict to store current-best IDs and scores for sp1
    sp2IDs = dict() #dict to store current-best IDs and scores for sp2
  

  # parse SAM (stdin), extract best alignment scores, store in dict
    print("splitByScore_v3 run started " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    print(vString)
    print("Chosen species is " + sp1 + ", alternate species is " + sp2)
  
  # parse SAM (stdin), extract best alignment scores, store in dict
    f = sys.stdin
    start = time.time()
    for line in f:
        if line[0] == '@':
            continue #skip header
        ReadCount += 1 #found a read line
        spl = line.rstrip().split('\t') #split into fields
        flag = int(spl[1])
        if flag & 0x4: 
            UnmappedCount += 1

    #read is mapped
        GoodCount += 1
        # retrieve alignment score
        aScore = getTag('AS', spl[11:])
        if not aScore:
            MissingAS += 1
            continue
        aScore = int(aScore)
    
        readGenome = spl[2][:len(sp1)]
        if readGenome == sp1:  #this read is from sp1
            if spl[0] in sp1IDs:  #seen this read in sp1 before
                if aScore > sp1IDs[spl[0]]: # this read has better AS
                    sp1IDs[spl[0]] = aScore  #replace with better AS
                    continue
                else:
          # previous read had equal or better AS -- do nothing
                    continue
            elif spl[0] in sp2IDs: #was is seen in sp2 before?
                if aScore > sp2IDs[spl[0]]:  #this new sp1 read is better than the old sp2 read
                    del sp2IDs[spl[0]]  #remove from sp2 dict
                    sp2Count += -1  #decrement sp2 count
                    sp1IDs[spl[0]] = aScore  #add to sp1 dict
                    sp1Count += 1  #increment sp1 count
                    continue
                else: # score in sp2 is better than in sp1, do nothing
                    continue
            else:
                #not seen this before - add it to sp1 dict
                sp1IDs[spl[0]] = aScore
                sp1Count += 1
        elif readGenome == sp2:  #this read is from sp2
            if spl[0] in sp2IDs:  #seen this read in sp2 before
                if aScore > sp2IDs[spl[0]]: #this read has better score
                    sp2IDs[spl[0]] = aScore  # replace with better AS
                    continue
                else:
                    # previous read had equal or better AS -- do nothing
                    continue
            elif spl[0] in sp1IDs:  # was this seen in sp1 before?
                if aScore > sp1IDs[spl[0]]: #this new sp2 read is better than the old sp1 read
                    del sp1IDs[spl[0]] #remove from sp1 dict
                    sp1Count += -1
                    sp2IDs[spl[0]] = aScore #add to sp2 dict
                    sp2Count += 1
                    continue
                else: # score in sp1 is better than in sp2 -- do nothing
                    continue
            else: #not seen this before, add it to sp2 dict
                sp2IDs[spl[0]] = aScore
                sp2Count += 1
  
    f.close()  #done with stdin
    
    end = time.time()
    print("Parsing SAM output took " + '{:.2f}'.format(end - start) + " seconds.")
      
      #Write status of ID matching to stdout or redirected to oFname
    print("Total " + numFmt(ReadCount) + " reads, " + numFmt(GoodCount) + " are aligned" + numFmt(UnmappedCount) + " are unmapped")
    print("number that miss AS", MissingAS)
    print(numFmt(sp1Count) + " reads matching " + sp1)
    print(numFmt(sp2Count) + " reads matching " + sp2)
    sys.stdout.flush() #to empty stdout buffer to file
      
      # Find paired-end fastq files.
    fastqs = [f for f in os.listdir(fastqDir) if f.endswith(fastqExt)]
    R1 = ''.join([f for f in fastqs if (readEnds[0] in f) & (SEL_STRING in f)])  
    R2 = ''.join([f for f in fastqs if (readEnds[1] in f) & (SEL_STRING in f)])
    fastqs = [ R1,R2 ]
    genDicts = [ sp1IDs, sp2IDs ]
    for genID,gen in enumerate(genomes):
      	for endID,readEnd in enumerate(readEnds):
      		fastqOut(gen,genDicts[genID],fastqDir,fastqs[endID])
    
    print("splitByScore_v3 run ended " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    
    sys.exit(0)
  
if __name__ == '__main__':
    main()
import sys,os
from Bio import SeqIO
import datetime
import time
import argparse

# Steps:
# 1. Read SAM input and parse into:
#    a. dictionary of read numbers pointing to tuple (immutable) required aligment info
#    b. list of unique read IDs pointing to list (updatable) of rows in dict
# 2. Report counts
# 3. Decide best genome:
#    a. If ID matches only one genome, that's it
#    b. If ID has multiple genome matches, choose primary alignment
#    c. (future) if any reads remaining, choose best alignment score
# 4. Output tables for splitting fastq files in a separate program 

# Global variables with standard naming convension
VERSION = '1.0'
DESCRIPTION = '''\
    The goal is to separate BULK RNAseq fastq files into two species specific subsets
    We iterate through the mixed-genome SAM file to and choose the best score by genome.
    Read IDs matching both genomes are counted.
	Open SAM input, parse into a dictionary of alignment info and
	a dictionary of read IDs (pointing to a tuple of alignment records).
	'''
def checkOptions():
    parser = argparse.ArgumentParser(usage='samtools view <BAM> | python byPrim.py',description=DESCRIPTION)
    parser.add_argument('-v','--version', action='version', version=VERSION)
    parser.add_argument('-g','--genome',dest="sp1",default='hg38',
                    help='Selected genome identifier (default hg38)')
    parser.add_argument('-a','--altgenome',dest="sp2",default='mm10',
                    help='Non-matching genome identifier (default mm10)')
    parser.add_argument('-s','--sample',dest="sample_name",required=True,default=None,
                    help='Sample ID for insertion into file name')

    args = parser.parse_args()
    #if args.sample_name:
        #args.sample_name = args.sample_name.upper()
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
	'''
	Format large numbers with commas for string output.
	'''
	return '{:,}'.format(value)

def unique(lis):
	'''
	Return list of unique values fron a list
	'''
	u = []
	for x in lis:
		if x not in u:
			u.append(x)
	return(u)

def bestScore(lis):
	'''
	Return highest alignment score from a list
	'''
	b = -999 #smaller than any score
	for x in lis:
		if x > b:
			b = x
	return(b)

def numPrimaries(lis):
	'''
	Return count of isPrimary == True
	'''
	c = 0
	for x in lis:
		if x:
			c += 1
	return(c)

def bestPrim(p,s):
	prims = []
	b = -999
	j = -1
	for i in range(len(p)): #iterates through indexes in p
		if p[i]:
			prims.append(i) #add index to list of primary alignments (should be 2 of them)
	for prim in prims:
		if s[prim] > b: #check both indices to see if they're the better score
			b = s[prim] #this has better score
			j = prim #this is the prims index to return
	return(j)
		

def main():
	
	opt = checkOptions()
	# all options are accessible as opt.var
	# print(opt)
	sp1 = opt.sp1
	sp2 = opt.sp2
	sampleName = opt.sample_name
  
	# initialize counts
	ReadCount = 0
	GoodCount = 0
	sp1Count = 0
	sp2Count = 0
	
	#Create empty dictionaries
	samTable = dict() # each entry will be GoodCount -> (genome,aScore,isPrim)
	readIDs = dict() # each entry will be readID -> [ rows in samTable ]
	
	#Start log output written to stdout
	print("byPrim started " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
	print("Chosen species is " + sp1 + ", alternate species is " + sp2)
	
	#Open sam file and parse into tables
	f = sys.stdin
	start = time.time()
	for line in f:
		if line[0] == '@':
			continue #skip header
		ReadCount += 1 #Found a read
		t = line.rstrip().split('\t') #split line into tokens
		flag = int(t[1])
		if flag & 0x4:
			continue #Read is unmapped
		GoodCount += 1 #This read is mapped
		genome = t[2][:t[2].index('_')] #slice chars before _ from third token
		aScore = getTag('AS',t[11:])
		if not aScore:
			print('missing AS: ' + t[0] + '\n')
			continue #check this -- any missing?
		aScore = int(aScore)
		if flag & 0x100:
			isPrim = False
		else:
			isPrim = True # True if this read is a primary read (not flag 256, 'not primary alignment')
		samTable[GoodCount]=(genome,aScore,isPrim) #only store needed fields: genome, aScore, isPrim
		if t[0] in readIDs.keys():  #we've seen this ID before
			readIDs[t[0]].append(GoodCount)
		else:
			readIDs[t[0]] = [GoodCount] #key points to set of row numbers as tuple
	
	f.close()
	end = time.time()
	print("Parsing SAM output took " + '{:.2f}'.format(end - start) + " seconds.")
	print("Found " + numFmt(ReadCount) + " records in SAM input.")
	print(numFmt(GoodCount) + " are mapped.")
	print("Stored " + numFmt(len(readIDs)) + " unique read IDs.")
	print("Pointing to " + numFmt(len(samTable)) + " alignment records.")
	print("\nParsing tables...")
	
	#build string for tabulated log
	logString = sampleName + '\t' + str(ReadCount) + '\t' + str(GoodCount) + '\t' + str(len(readIDs)) + '\t' + str(len(samTable)) + '\t'
	
	#Loop through read IDs
	#Count sam records and number of unique genomes
	#if unique genomes == one, choose this genome and remove these records from dicts
	start = time.time()
	
	#set counters
	NumGen1 = 0
	NumGen2 = 0
	NumMixGen1 = 0
	NumMixGen2 = 0
	
	#create lists of records to be removed
	readIDs2del = []
	samTable2del = []
	
	f = open(sampleName + "_table.tsv",'w') #output file to store decision table
	for rID in readIDs: #loop through keys of readIDs dictionary
		nSamRecs = len(readIDs[rID]) #number of elements in value list
		theseSamRecs = [ samTable[x] for x in readIDs[rID] ] #returns list of lists for rID
		genomes = [ tup[0] for tup in theseSamRecs ]
		scores = [ tup[1] for tup in theseSamRecs ]
		primaries = [ tup[2] for tup in theseSamRecs ]
		
		#time for decisions:
		#if read points to one unique genome, we're done, code is uniq for only one genome
		if len(unique(genomes)) == 1:
			#write output to output table file then remove records from dicts
			if unique(genomes)[0] == sp1:
				NumGen1 += 1
			elif unique(genomes)[0] == sp2:
				NumGen2 += 1
			f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(rID,'uniq',len(readIDs[rID]),unique(genomes)[0],bestScore(scores),numPrimaries(primaries)))
			GoodCount -= nSamRecs #decrement count by the number of alignments selected by rID
			for key in readIDs[rID]: 
				samTable2del.append(key)
			readIDs2del.append(rID)
		else:
			#more than one genome for this read ID, code is prim for choice based on primary alignment
			whichIsPrim = bestPrim(primaries,scores) #returns index of primary with best score
			if scores[whichIsPrim] < bestScore(scores):
				print("Found primary alignment with score not best: " + rID + " with primary score=" + str(scores[whichIsPrim]) + " but best=" + str(bestScore(scores)))
			if genomes[whichIsPrim] == sp1:
				NumMixGen1 += 1
			elif genomes[whichIsPrim] == sp2:
				NumMixGen2 += 1
			f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(rID,'prim',len(readIDs[rID]),genomes[whichIsPrim],bestScore(scores),numPrimaries(primaries)))
			GoodCount -= nSamRecs #decrement count by the number of alignments selected by rID
			for key in readIDs[rID]: 
				samTable2del.append(key)
			readIDs2del.append(rID)
	f.close()
	
	print("After storing decision output in " + sampleName + "_table.tsv")
	print("Remaining readIDs: " + numFmt(GoodCount))
	print("Number uniquely matching " + sp1 +": " + numFmt(NumGen1))
	print("Number uniquely matching " + sp2 +": " + numFmt(NumGen2))
	print("Total unique matches: " + numFmt(NumGen1 + NumGen2))
	print("Number mixed: " + numFmt(NumMixGen1 + NumMixGen2))
	print("Of mixed: " + numFmt(NumMixGen1) + " match " + sp1 + " and " + numFmt(NumMixGen2) + " match " + sp2)
	print("Total ouput: " + numFmt(NumGen1 + NumGen2 + NumMixGen1 + NumMixGen2))
	end = time.time()
	print("Deciding tables and writing output took " + '{:.2f}'.format(end - start) + " seconds.")
	
	for rID in readIDs2del:
		del readIDs[rID]
	for key in samTable2del:
		del samTable[key]
	
	#diagnostics if any reads end up here
	f = open(sampleName + '_remaining_dict.txt','wt')
	f.write(str(samTable))
	f.close()
	
	#final step - write a table-format log for concatenation with other samples
	f = open(sampleName + "_log.tab",'wt')
	f.write(logString + str(NumGen1) + '\t' + str(NumGen2) + '\t' + str(NumMixGen1) + '\t' + str(NumMixGen2) + '\t' + str(GoodCount) + '\n')
	f.close()

	
	sys.exit(0)
	
if __name__ == '__main__':
	main()

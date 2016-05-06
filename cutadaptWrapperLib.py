import os,sys
from cPickle import load,dump
import re

def trashExtractor(inp):
        """ retrieve Overrepresented sequences from a fastqc_data file (i.e. polyA and adapters) """
        with open(inp) as f:
                flag = 0
                for line in f:
                        line = line.strip()
                        if not flag:
                                if not line.startswith('>>Overrepresented sequences'): continue
                                flag = 1
                        else:
                                if line.startswith('>>END_'): return
                                if not line.startswith('#'): yield line


def getAdaptersAndPolyAs_old(trash,adaptersDict=None):
        """ return list of sequence for makeCutAdaptCL from fastqc_data file (old vrsn) """
	if not isinstance(adaptersDict,dict):
		try: adaptersDict=load(open('adaptersIlluminaTruSeq.pkl'))
		except: raise Exception("can't find adaptersIlluminaTruSeq.pkl in your local PATH")
        adapters,polyAs,universal = set(),[],None
	m = re.compile("TruSeq Adapter, Index \d+")
	for line in trash:
		match = m.search(line)
		if match: adapters.add(adaptersDict[match.group().lower()])
		elif "Illumina Single End PCR Primer 1" in line: universal="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
		else: polyAs.append(line.split('\t')[0])
        return adapters,polyAs,universal


def getAdaptersAndPolyAs(trash,adaptersDict=None):
        """ return list of sequence for makeCutAdaptCL from fastqc_data file """
	if not isinstance(adaptersDict,dict):
		try: adaptersDict=load(open('adaptersIlluminaTruSeq.pkl'))
		except: raise Exception("can't find adaptersIlluminaTruSeq.pkl in your local PATH")
        adapters,polyAs,universal = set(),[],None
	m = re.compile("TruSeq Adapter, Index \d+")
	for line in trash:
		match = m.search(line)
		if match: adapters.add(line.strip().split()[0])
		elif "Illumina Single End PCR Primer 1" in line: universal="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
		else: polyAs.append(line.split('\t')[0])
        return adapters,polyAs,universal

def makeAdaptersPkl(inp="./adapters.fasta"):
	d = {r.description:str(r.seq) for r in parse(inp,'fasta')}
	with open('adaptersIlluminaTruSeq.pkl', 'wb') as fp: dump(d,fp)
	return

def makeCutAdaptCL(inp1,inp2,adapters,universal=None):
        """ prepare Commmand line for cutadapt. Retrieve primers using the fastqc_data files and the adapters.fasta file 
        Default cutadapt call is (for TrueSeq; borrowed from cutadapt manual):
                cutadapt \
                    -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
                    -o trimmed.1.fastq.gz -p trimmed.2.fastq.gz \
                    reads.1.fastq.gz reads.2.fastq.gz
        output: a string with the command line for the cutadapt call
                """
        out1,out2 = os.path.join(os.path.dirname(inp1),'trimmed_'+os.path.basename(inp1)),os.path.join(os.path.dirname(inp2),'trimmed_'+os.path.basename(inp2))
        # if universal==None: universal="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
        adapts = " ".join(["-a %s" %i for i in adapters])
        if universal: commLine = "cutadapt %s -A %s -o %s -p %s %s %s" %(adapts,universal,out1,out2,inp1,inp2)
	else: commLine = "cutadapt %s -o %s -p %s %s %s" %(adapts,out1,out2,inp1,inp2)
        return commLine,out1,out2


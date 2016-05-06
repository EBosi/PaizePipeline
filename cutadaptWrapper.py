from pythonFxns import *
from cPickle import load,dump
import sys,os

if __name__ == "__main__":
	__usage__ = " python cutadaptWrapper.py file1 file2 ---> call line for cutadapt "
	inps = sys.argv[1:]
	if not 'adaptersIlluminaTruSeq.pkl' in os.listdir('.'): makeAdaptersPkl()
	f1,f2 = inps
	extractedFromF1,extractedFromF2 = getAdaptersAndPolyAs(trashExtractor(f1)),getAdaptersAndPolyAs(trashExtractor(f2))
	adapters,polyAs,universal = set.union(extractedFromF1[0],extractedFromF2[0]),set(extractedFromF1[1] + extractedFromF2[1]),extractedFromF1[-1] or extractedFromF2[-1]
	adapters = list(adapters) + list(polyAs)
	print makeCutAdaptCL(inp1,inp2,universal=universal,*adapters)

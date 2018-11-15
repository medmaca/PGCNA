#!/usr/bin/env python
"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Matthew Care

A memory efficient implementation of PGCNA, that requires very little RAM
to process very large expression data-sets.
"""
from __future__ import print_function

import sys
import os
import argparse
import shutil
from datetime import datetime
from collections import defaultdict
import re
import glob
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import h5py


##----------------------------------------Create argparse argments----------------------------------------##
parser = argparse.ArgumentParser()
#  Required
parser.add_argument("-w","--workFolder",help="Work Folder [REQUIRED]",default=None)
parser.add_argument("-d","--dataF",help="Expression data file path [REQUIRED]",default=None)

#  Optional
parser.add_argument("-s","--fileSep",help="Separator used in expression file [\t]",default="\t")
parser.add_argument("--headerL",help="Number of header lines before expression values [1]",default=1,type=int)
parser.add_argument("-f","--retainF",help="Retain gene fraction -- keeping most variant genes [0.8]",type=float,default=0.8)
parser.add_argument("-e","--edgePG",help="Edges to keep per gene [3] -- Highly recommend leaving as default",type=int,default=3)

#  Not often changed
parser.add_argument("--outF",help="Root output folder [PGCNA]",default="PGCNA")
parser.add_argument("--corrMatF",help="Correlation matrix folder [CORR_MATRIX]",default="CORR_MATRIX")
parser.add_argument("--gephiF",help="Folder to store files for Gephi [GEPHI]",default="GEPHI")
parser.add_argument("--fastUF",help="Fast Unfolding folder [FAST_UNFOLD]",default="FAST_UNFOLD")

#  Flags
parser.add_argument("--noFastUF",help="Don't try to run Fast Unfolding Clustering, but complete everything else [False] -- Flag",action="store_true")
parser.add_argument("--usePearson",help="Use Pearson Correlation instead of Spearman [False] -- Flag",action="store_true")
parser.add_argument("--keepBigF",help="Keep big HDF5 files after finishing [False] -- Flag",action="store_true")

#  Correlation related
parser.add_argument("--corrChunk",dest="corrChunk",help="Size of chunk (rows) to split correlation problem over [5000] -- Higher will speed up correlation calculation at cost of RAM",default=5000,type=float)

#  FastUnfolding specific
parser.add_argument("-n","--fuNumber",dest="fuRunNum",help="Number of times to run [100]",default=100,type=float)
parser.add_argument("-r","--fuRetain",dest="fuRetainNum",help="Retain top [1] clusterings",default=1,type=float)
parser.add_argument("--fuRenumberStart",dest="fuRenumberStart",help="FU clustering re-numbering start [1]",default=1,type=int)
parser.add_argument("--tOutFold",dest="tOutFold",help="Trees Ouput folder [Trees]",default="Trees")
parser.add_argument("--fTOutFold",dest="fTOutFold",help="Final retained Trees Ouput folder [TreesF]",default="TreesF")
parser.add_argument("--cOutFold",dest="cOutFold",help="Clusters out folder [Clusters]",default="Clusters")
parser.add_argument("--cOutTxtFold",dest="cOutTxtFold",help="Clusters out folder - mapped back to genes [ClustersTxt]",default="ClustersTxt")
parser.add_argument("--cOutListsFold",dest="cOutListsFold",help="Clusters out folder - split into lists [ClustersLists]",default="ClustersLists")


###########################################################################################
args = parser.parse_args()
if not args.workFolder or not args.dataF:
	print("\n\nNeed to specifiy REQUIRED variables see help (-h)")
	sys.exit()

###########################################################################################
OUT_FOLDER = os.path.join(args.workFolder,args.outF,"EPG" + str(args.edgePG))
if not os.path.exists(OUT_FOLDER):
	os.makedirs(OUT_FOLDER)
###########################################################################################
###########################################################################################
if not args.noFastUF:
	def testCommandOut(outR,testFail):
		if re.search("command not found",outR):  # /bin/sh: louvainn: command not found
			return True
		elif re.search("not recognized as an internal or external command",outR):
			return True

	#  Make sure that if Fast Unfolding is to be run it exists on the path
	testCommand = "louvain -h"
	p = Popen(testCommand,shell=True,stdout=PIPE,stderr=STDOUT)
	stdout,stderr = p.communicate()

	testCommand = "convert"
	p = Popen(testCommand,shell=True,stdout=PIPE,stderr=STDOUT)
	stdout2,stderr = p.communicate()

	testCommand = "hierarchy"
	p = Popen(testCommand,shell=True,stdout=PIPE,stderr=STDOUT)
	stdout3,stderr = p.communicate()

	testFail = False
	testFail = testCommandOut(stdout,testFail)
	testFail = testCommandOut(stdout2,testFail)
	testFail = testCommandOut(stdout3,testFail)

	if testFail:
		print("\nCouldn't run Fast unfold commands (louvainn, convert and hierarchy) must be in your path/environment (to skip clustering step use --noFastUF flag), see readme for more information")
		sys.exit()
###########################################################################################
#  Tidy input
#  Decode fileSep for passed "\t"
args.fileSep = args.fileSep.decode('string-escape')
args.corrChunk = int(round(args.corrChunk,0))
args.fuRunNum = int(round(args.fuRunNum,0))
args.fuRetainNum = int(round(args.fuRetainNum,0))
##----------------------------------------LOGGER----------------------------------------##


class multicaster(object):
	def __init__(self, filelist):
		self.filelist = filelist

	def write(self, str):
		for f in self.filelist:
			f.write(str)


def concatenateAsString(joinWith,*args):
	temp = [str(x) for x in args]
	return joinWith.join(temp)

print("##----------------------------------------",str(datetime.now()),"----------------------------------------##",sep="")
#  print out settings
settingInf = concatenateAsString(
	"\n",
	"##----------------------------------------Arguments----------------------------------------##",
	"#  Required",
	"WorkFolder [-w,--workFolder] = " + args.workFolder,
	"Expresion data file path [-d,--dataF] = " + args.dataF,
	"\n#  Optional",
	"Separator used in expression file [-s, --fileSep] = " + args.fileSep,
	"Number of header lines in expression file [--headerL] = " + str(args.headerL),
	"Retain gene fraction [-f, --retainF] = " + str(args.retainF),
	"Edges to retain per gene [-e, --edgePG] = " + str(args.edgePG),
	"\n#  Not changed often"
	"Root output folder [--outF] = " + args.outF,
	"Correlation matrix folder [--corrMatF] = " + args.corrMatF,
	"Gephi files folder [--gephiF] = " + args.gephiF,
	"Fast unfolding folder [--fastUF] = " + args.fastUF,
	"\n#  Flags",
	"Don't run Fast Unfolding clustering [--noFastUF] = " + str(args.noFastUF),
	"Use Pearson Correlation [--usePearson] = " + str(args.usePearson),
	"Keep big HDF5 files after run [--keepBigF] = " + str(args.keepBigF),
	"\n#  Correlation Options",
	"Chunk size (rows) to split correlation over [--corrChunk]:\t" + str(args.corrChunk),
	"\n#  Fast-Unfolding Specific",
	"Cluster output folder -- Clusters split into lists [--cOutListsFold]:\t" + str(args.cOutListsFold),
	"Run number [-n, --fuNumber]:\t" + str(args.fuRunNum),
	"Retain top number of clusterings [-r, --fuRetain]:\t" + str(args.fuRetainNum),
	"Fast-Unfolding clustering renumber start [--fuRenumberStart]:\t" + str(args.fuRenumberStart),
	"Tree output folder [--tOutFold]:\t" + str(args.tOutFold),
	"Final retained tree output folder [--fTOutFold]:\t" + str(args.fTOutFold),
	"Cluster output folder [--cOutFold]:\t" + str(args.cOutFold),
	"Cluster output folder -- mapped back to genes [--cOutTxtFold]:\t" + str(args.cOutTxtFold),
)

settingsF = open(os.path.join(OUT_FOLDER,str(datetime.now()).replace(" ","-").replace(":",".")[:-7] + "_CorrelationSettingsInfo.txt"),"w")
print(settingInf,file=settingsF)
settingsF.close()
###########################################################################################


##----------------------------------------Methods----------------------------------------##

def make_key_naturalSort():
	"""
	A factory function: creates a key function to use in sort.
	Sort data naturally
	"""
	def nSort(s):
		convert = lambda text: int(text) if text.isdigit() else text
		alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

		return alphanum_key(s)
	return nSort


def generateCorrMatrix(workF,dataF,retainF=0.8,corrMatF="CORR_MATRIX",usePearson=False,corrChunk=5000,splitBy="\t",headerL=1,decimalP=3,decimalPpVal=5):
	"""
	Memory efficient method for generating all pairwise correlations of genes (rows) across a set of samples (columns).  Uses HDF5 files to greatly
	reduce memory useage, keeping most data residing on the hard disk.
	"""
	print("\n\nGenerate correlations for expression data:")



	# Create subfolder
	outF = os.path.join(workF,corrMatF)
	if not os.path.exists(outF):
		os.makedirs(outF)
	
	#############################
	# Read in data file
	seenGenes = {}
	genes = []
	expression = []
	header = headerL
	fileN = os.path.basename(dataF)
	print("\tReading expression data file:",fileN)

	for i, line in enumerate(open(dataF)):
		cols = line.rstrip().split(splitBy)

		if header:
			header -= 1
			continue

		genes.append(cols[0])
		if not cols[0] in seenGenes:
			seenGenes[cols[0]] = 1
		else:
			print("\n\nERROR: Duplicate gene (",cols[0],") in expression data (",fileN,"), expects unique identifiers!\nPlease remove duplicates and re-run\n",sep="")
			sys.exit()
		expression.append(cols[1:])

	print("\t\tTotal genes = ",len(genes))
	#############################

	#############################
	print("\tGenerate Numpy matrices")
	#  Move genes to numpy array, so can be sorted along with the rest
	genesNP = np.empty((len(genes)),dtype=str)
	genesNP = np.array(genes)

	#  Create numpy array
	nrow = len(expression[0])
	ncol = len(expression)
	npA = np.empty((ncol,nrow),dtype=object)
	for i,eSet in enumerate(expression):
		npA[i] = eSet

	#  Calculate SD for each gene
	npSD = np.std(npA.astype(np.float64),axis=1)

	#  Sort by SD
	sortI = np.argsort(npSD)[::-1]  # Sort ascending
	#  Sort matrix by index
	npA = npA[sortI,:]

	#  Sort genes by index
	genesNP = genesNP[sortI]

	#  Only retain top fract #
	#  Work out index to cut at
	cutPos = int(round(retainF * genesNP.shape[0],0))
	print("\t\tRetain Fraction (",retainF,") Cut at: ",cutPos,sep="")
	print("\t\t\tPre-cut shape:",npA.shape)

	#  Retain most variant
	npA = npA[:cutPos,].astype(np.float64)
	genesNP = genesNP[:cutPos]
	print("\t\t\tPost-cut shape:",npA.shape)
	
	#  Output genes
	outMatrixN = os.path.join(outF,fileN + "_RetainF" + str(retainF) + "_Genes.txt")
	np.savetxt(outMatrixN,genesNP,fmt="%s",delimiter="\t")
	
	#######################################
	# Calculate correlations
	print("\tCalculating correlations using HDF5 file to save memory")
	rowN, colN = npA.shape
	
	outMatrixHDF5 = os.path.join(outF,fileN + "_RetainF" + str(retainF) + ".h5")
	with h5py.File(outMatrixHDF5, "w") as f:
		h5 = f.create_dataset("corr", (rowN,rowN), dtype="f8")

	# # Load into memory
	print("\t\tCreating HDF5 file")
	h5 = h5py.File(outMatrixHDF5,'r+')
	
	if not usePearson:
		# Note this isn't a perfect Spearman, as ties are not dealt with, but it is fast, and given we're only
		# Retaining the most correlated edges will have almost no effect on the output.
		npA = npA.argsort(axis=1).argsort(axis=1).astype(np.float64)
	
	# subtract means from the input data
	npA -= np.mean(npA, axis=1)[:,None]

	# normalize the data
	npA /= np.sqrt(np.sum(npA*npA, axis=1))[:,None]

	# Calculate correlations per chunk
	print("\t\tCalculating correlations for Chunk:")
	for r in range(0, rowN, corrChunk):
		print("\t\t\t",r)
		for c in range(0, rowN, corrChunk):
			r1 = r + corrChunk
			c1 = c + corrChunk
			chunk1 = npA[r:r1]
			chunk2 = npA[c:c1]
			h5["corr"][r:r1, c:c1] = np.dot(chunk1, chunk2.T)

	# # Write output out for debugging
	# finalCorr = np.copy(h5["corr"][:])
	# outMatrixN = os.path.join(outF,fileN + "_RetainF" + str(retainF) + "_Corr")
	# stringFormat = "%." + str(decimalP) + "f"
	# # stringFormat2 = "%." + str(decimalPpVal) + "f"
	# np.savetxt(outMatrixN + ".txt",finalCorr,fmt=stringFormat,delimiter="\t")

	#  Return path to HDF5 file and genes matrix
	return outMatrixHDF5, genesNP


def reduceEdges(workF,dataF,gephiF,corrh5,genesM,retainF=0.8,edgePG=3,printEveryDiv=10,corrMatF=None,keepBigF=False):
	"""
	Reduce edges in correlation matrix, only retaining edgePG maximum correlated genes per row
	"""
	fileN = os.path.basename(dataF)

	print("\tLoad HDF5 file")
	#  Load old HDF5
	h5 = h5py.File(corrh5,'r+')
	rowN, colN = h5["corr"].shape

	printEvery = int(round(rowN/float(printEveryDiv),0))
	print("\nReduces edges to (",edgePG,") per gene:",sep="")
	printE = printEvery
	tell = printE
	count = 0

	print("\tWorking (report every 1/",printEveryDiv,"th of total):",sep="")
	for i, row in enumerate(h5["corr"]):
		#  Inform user of position
		count += 1
		if count == printE:
			printE = printE + tell
			print("\t\t",count)

		sortI = np.argsort(row)[::-1]  # Descending
		
		# Need to find last index with same value (in case have a run of exactly the same correlation)
		tmp = row[sortI]
		maxI = edgePG  # Set off from last but one element
		lastElement = tmp[maxI]
		try:
			while True:
				maxI += 1
				if maxI == row.shape[0]:
					break
				if not tmp[maxI] == lastElement:
					maxI -= 1
					break
		except IndexError:
			print("Error:",row.shape,tmp)
			sys.exit()

		# Set all other values to zero
		t = h5["corr"][i]  # Need to create a temporary copy, as h5 doesn't like writing to a sorted range.
		t[sortI[maxI+1:]] = 0
		h5["corr"][i] = t

	# Write output out for debugging
	# finalCorr = np.copy(h5EPG["corr"][:])
	# outMatrixN = os.path.join(outF,fileN + "_RetainF" + str(retainF) + "_NonSym")
	# stringFormat = "%." + str(3) + "f"
	# # stringFormat2 = "%." + str(decimalPpVal) + "f"
	# np.savetxt(outMatrixN + ".txt",finalCorr,fmt=stringFormat,delimiter="\t")

	###########################################################################################
	
	print("\nGenerating files for Gephi network visualization tool")
	#  Create subfolder
	outF = os.path.join(workF,gephiF)
	if not os.path.exists(outF):
		os.makedirs(outF)

	#  First output list of genes (nodes)
	print("\tCreate node file")
	nodesFP = os.path.join(outF,fileN + "_RetainF" + str(retainF) + "_EPG" + str(edgePG) + "_Nodes.tsv")
	nodesFile = open(nodesFP,"w")
	print("Id\tLabel\tCluster",file=nodesFile)
	for gene in genesM:
			print(gene,gene,"NA",file=nodesFile,sep="\t")
	nodesFile.close()

	#  Second, output list of edges
	print("\tCreate edges file (report every 1/",printEveryDiv,"th of total):",sep="")
	edgesFP = os.path.join(outF,fileN + "_RetainF" + str(retainF) + "_EPG" + str(edgePG) + "_Edges.tsv")
	edgesFile = open(edgesFP,"w")
	print("Source\tTarget\tWeight\tType\tfromAltName\ttoAltName",file=edgesFile)
	
	# Finally output edges
	printE = printEvery
	tell = printE
	count = 0
	seen = defaultdict()
	printedLines = 0

	for i,row in enumerate(h5["corr"]):
		for j,item in enumerate(row):
			if not (i == j):  # Don't output self edges...
				geneO = genesM[i]
				geneT = genesM[j]
				if not (geneO + "-@@@-" + geneT) in seen:
					#  Output info
					if not item == 0:
						print(geneO,geneT,item,"undirected",geneO,geneT,sep="\t",file=edgesFile)
						printedLines += 1

						#  Store the fact that we've see this and it's equivalent pair
						seen[geneO + "-@@@-" + geneT] = 1
						seen[geneT + "-@@@-" + geneO] = 1

		count += 1
		if count == printEvery:
			printEvery = printEvery + tell
			print("\t\t",count,"Genes",", Edges:",printedLines)

	edgesFile.close()
	print("\n\t\tTotal printed (",printedLines,") edges",sep="")

	if not keepBigF:
		h5.close()
		shutil.rmtree(os.path.join(workF,corrMatF))

	return edgesFP, nodesFP


def prepareFastUnfold(workF,dataF,fastUF,edgesFP,retainF=0.8,edgePG=3,splitBy="\t",header=1):
	print("\nPreparing files for clustering using Fast Unfolding method")

	fileN = os.path.basename(dataF)
	#  Create subfolder
	outF = os.path.join(workF,fastUF,fileN + "_retainF" + str(retainF) + "_EPG" + str(edgePG))
	if not os.path.exists(outF):
		os.makedirs(outF)

	nodePairs = defaultdict(dict)
	naturalKey = make_key_naturalSort()

	nameToInt = {}
	nameIntCount = 0

	edgeNum = 0
	for line in open(edgesFP):
		cols = line.rstrip().split(splitBy)

		if header:
			header -= 1
			continue

		gene1 = cols[0]
		gene2 = cols[1]
		edgeNum += 1

		#  Convert from name to int
		if gene1 in nameToInt:
			gene1 = nameToInt[gene1]
		else:
			nameToInt[gene1] = str(nameIntCount)
			gene1 = str(nameIntCount)
			nameIntCount += 1

		if gene2 in nameToInt:
			gene2 = nameToInt[gene2]
		else:
			nameToInt[gene2] = str(nameIntCount)
			gene2 = str(nameIntCount)
			nameIntCount += 1
		
		#################################

		nodePairs[gene1][gene2] = cols[2]

	print("\tTotal edges =",edgeNum)

	outFilePath = os.path.join(outF,"sym_edges.txt")
	outFile = open(outFilePath,"w")
	for node in sorted(nodePairs,key=naturalKey):
		for node2 in sorted(nodePairs[node],key=naturalKey):
			print(node,node2,nodePairs[node][node2],sep="\t",file=outFile)
	outFile.close()

	#  Write out mapping between genes and ints
	outFilePath2 = os.path.join(outF,"GeneIntMap.txt")
	outFile2 = open(outFilePath2,"w")
	for gene in sorted(nameToInt,key=naturalKey):
		print(gene,nameToInt[gene],sep="\t",file=outFile2)
	outFile2.close()

	return outFilePath, outFilePath2


def runFastUF(symEdgesF,geneIntMapF,tOutFold,fTOutFold,cOutFold,cOutTxtFold,cOutListsFold,runNum=10,retainNum=1,printEveryDiv=10):
	MOD_SCORES = "modScores.txt"

	print("\nRunning Fast Unfolding")
	baseFold,symEdgesFile = os.path.split(symEdgesF)

	# Make folders
	T_OUT_FOLDER = os.path.join(baseFold,tOutFold)
	TF_OUT_FOLDER = os.path.join(baseFold,fTOutFold)
	C_OUT_FOLDER = os.path.join(baseFold,cOutFold)
	C_OUT_TXT_FOLDER = os.path.join(baseFold,cOutTxtFold)
	C_OUT_LISTS_FOLDER = os.path.join(baseFold,cOutListsFold)

	# Remove old
	if os.path.exists(TF_OUT_FOLDER):
		shutil.rmtree(TF_OUT_FOLDER)
	if os.path.exists(T_OUT_FOLDER):
		shutil.rmtree(T_OUT_FOLDER)
	if os.path.exists(C_OUT_FOLDER):
		shutil.rmtree(C_OUT_FOLDER)
	if os.path.exists(C_OUT_TXT_FOLDER):
		shutil.rmtree(C_OUT_TXT_FOLDER)
	if os.path.exists(C_OUT_LISTS_FOLDER):
		shutil.rmtree(C_OUT_LISTS_FOLDER)

	# Replace
	if not os.path.exists(T_OUT_FOLDER):
		os.makedirs(T_OUT_FOLDER)
	if not os.path.exists(TF_OUT_FOLDER):
		os.makedirs(TF_OUT_FOLDER)
	if not os.path.exists(C_OUT_FOLDER):
		os.makedirs(C_OUT_FOLDER)
	if not os.path.exists(C_OUT_TXT_FOLDER):
		os.makedirs(C_OUT_TXT_FOLDER)
	if not os.path.exists(C_OUT_LISTS_FOLDER):
		os.makedirs(C_OUT_LISTS_FOLDER)

	#####################################	
	
	#  Fast Unfold commands
	def convertBinary(iF):
		command = "convert -i " + iF + " -o " + iF[:-3]+"bin " + " -w " + iF[:-3] + "weights"
		p = Popen(command,shell=True,stdout=PIPE,stderr=STDOUT)
		stdout,stderr = p.communicate()

	def cluster(wF,iF,outFolder,runNum,PRINT_EVERY):
		modScores = []

		printEvery = PRINT_EVERY
		tell = printEvery
		count = 0

		for i in range(runNum):
			command = "louvain " + iF[:-3]+"bin -w " + iF[:-3] + "weights -l -1 > " + os.path.join(outFolder,str(i))
			p = Popen(command,shell=True,stdout=PIPE,stderr=STDOUT)
			stdout,stderr = p.communicate()
			modScores.append(stdout.rstrip())

			#  Notify position
			count += 1
			if count == printEvery:
				printEvery = printEvery + tell
				print("\t",count)

		outFile = open(os.path.join(wF,MOD_SCORES),"w")
		for i,score in enumerate(modScores):
			print(i,score,sep="\t",file=outFile)
		outFile.close()

	def selectClusterings(wF,modS,retainNumber,treeFold,finalTreeFold,splitBy="\t"):
		clusters = []
		scores = []

		for line in open(os.path.join(wF,modS)):
			cols = line.rstrip().split(splitBy)

			clusters.append(cols[0])
			scores.append(float(cols[1]))

		j = zip(scores,clusters)
		jSorted = sorted(j,reverse=True)
		
		scoresS,clustersS = zip(*jSorted)
		retained = clustersS[:retainNumber]

		for clst in retained:
			shutil.copy(os.path.join(treeFold,clst),os.path.join(finalTreeFold,clst))

	def outputClusters(finalTreeFold,clustFold):
		for path in glob.glob(os.path.join(finalTreeFold,"*")):
			fileName = os.path.basename(path)
			command = "hierarchy " + path
			p = Popen(command,shell=True,stdout=PIPE,stderr=STDOUT)
			stdout,stderr = p.communicate()

			bottomLevel = int(stdout.split("\n")[0].split(":")[-1].rstrip()) - 1
			
			command = "hierarchy " + path + " -l " + str(bottomLevel) + " > " + os.path.join(clustFold,fileName)
			p = Popen(command,shell=True,stdout=PIPE,stderr=STDOUT)
			stdout,stderr = p.communicate()

	def readMappings(mF,splitBy="\t"):

		intToGene = {}

		for line in open(mF):
			cols = line.rstrip().split(splitBy)

			if cols[0] in intToGene:
				print("ALREADY SEEN",cols[0],"SHOULD BE NON-REDUNDANT!")
				sys.exit()

			intToGene[cols[1]] = cols[0]

		return intToGene

	def processClusters(mapping,inFold,outFold,outListsFold,renumber_start=1,splitBy=" "):

		natKey = make_key_naturalSort()

		fileNumber = renumber_start
		for path in sorted(glob.glob(os.path.join(inFold,"*")),key=natKey):
			outFile = open(os.path.join(outFold,str(fileNumber) + ".csv"),"w")
			
			#  Create folder to store split gene lists in
			listFold = os.path.join(outListsFold,"Clust" + str(fileNumber))
			if not os.path.exists(listFold):
				os.makedirs(listFold)

			fileNumber += 1

			print("Id,Modularity Class",file=outFile)
			byClust = defaultdict(list)

			for line in open(path):
				cols = line.rstrip().split(splitBy)
				
				gene, clust = cols
				clust = str(int(clust) + 1)  # Start numbering from 1 not 0
				if gene in mapping:
					byClust[clust].append(mapping[gene])
				else:
					print("Gene (",gene,") is missing from mapping file!")

			for clst in sorted(byClust,key=natKey):
				tempLF = open(os.path.join(listFold,"M" + str(clst) + ".txt"),"w")
				for gene in sorted(byClust[clst],key=natKey):
					print(gene,clst,sep=",",file=outFile)
					print(gene,file=tempLF)
				tempLF.close()
			outFile.close()

	#####################################
	#####################################

	print("\tConvert file to binary")
	convertBinary(symEdgesF)

	print("\tRunning clusterings")
	printEvery = int(round(runNum/float(printEveryDiv),0))

	cluster(baseFold,symEdgesF,T_OUT_FOLDER,runNum,printEvery)

	print("\tRetrieving best",args.fuRetainNum,"clusterings")
	selectClusterings(baseFold,MOD_SCORES,retainNum,T_OUT_FOLDER,TF_OUT_FOLDER)

	print("\tOutput clusterings from bottom level")
	outputClusters(TF_OUT_FOLDER,C_OUT_FOLDER)

	#  Post-process results back to human readable files
	print("\tMap nodes back to genes and split into lists")
	mappings = readMappings(geneIntMapF)
	processClusters(mappings,C_OUT_FOLDER,C_OUT_TXT_FOLDER,C_OUT_LISTS_FOLDER,renumber_start=args.fuRenumberStart)


###########################################################################################

def main(finishT="Finished!"):

	#  Calculate correlations
	corrh5,genesM = generateCorrMatrix(OUT_FOLDER,args.dataF,retainF=args.retainF,corrMatF=args.corrMatF,usePearson=args.usePearson,corrChunk=args.corrChunk,splitBy=args.fileSep,headerL=args.headerL)
	
	#  Reduce edges
	edgesFP,nodesFP = reduceEdges(OUT_FOLDER,args.dataF,args.gephiF,corrh5,genesM,retainF=args.retainF,edgePG=args.edgePG,corrMatF=args.corrMatF,keepBigF=args.keepBigF)
	
	#  Prepare and run Fast Unfolding Of Communities in Large Networks (https://sourceforge.net/projects/louvain/files/louvain-generic.tar.gz/download)
	outFP1,outFP2 = prepareFastUnfold(OUT_FOLDER,args.dataF,args.fastUF,edgesFP,retainF=args.retainF,edgePG=args.edgePG)
	
	#  Run Fast Unfolding
	if not args.noFastUF:
		runFastUF(outFP1,outFP2,args.tOutFold,args.fTOutFold,args.cOutFold,args.cOutTxtFold,args.cOutListsFold,runNum=args.fuRunNum,retainNum=args.fuRetainNum)
	else:
		print("\nSkipping Fast-Unfolding clustering and downstream output as --noFastUF is set")
	print(finishT)


if __name__ == '__main__':
	main()
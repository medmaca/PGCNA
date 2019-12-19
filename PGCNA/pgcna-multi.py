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

This version is for processing multiple data-sets, via calculating a
median correlation matrix.
"""
from __future__ import print_function

import sys
import os
import argparse
import shutil
import string
from datetime import datetime
from collections import defaultdict
import re
import gzip
import glob
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import numpy.ma as ma
import scipy.stats as stats
import h5py


##----------------------------------------Create argparse argments----------------------------------------##
parser = argparse.ArgumentParser()
#  Required
parser.add_argument("-w","--workFolder",help="Work Folder [REQUIRED]",default=None)
parser.add_argument("-d","--dataF",help="Expression data folder path [REQUIRED]",default=None)

#  Optional
parser.add_argument("-m","--metaInf",help="File containing information about data files (must be in same folder as --dataF) [#FileInfo.txt]",default="#FileInfo.txt")
parser.add_argument("-s","--fileSep",help="Separator used in expression file(s) [\t]",default="\t")
parser.add_argument("-f","--retainF",help="Retain gene fraction -- keeping most variant genes [0.8]",type=float,default=0.8)
parser.add_argument("-g","--geneFrac",help="Fraction of files a gene needs to be present in to be included in median corr matrix [1/3]",type=float,default=1/3.0)
parser.add_argument("-e","--edgePG",help="Edges to keep per gene [3] -- Highly recommend leaving as default",type=int,default=3)

#  Not often changed
parser.add_argument("--outF",help="Root output folder [PGCNA]",default="PGCNA")
parser.add_argument("--corrMatF",help="Correlation matrix folder [CORR_MATRIX]",default="CORR_MATRIX")
parser.add_argument("--corrMatFS",help="Folder for single gene correlation files [CORR_MATRIX_SG]",default="CORR_MATRIX_SG")
parser.add_argument("--gephiF",help="Folder to store files for Gephi [GEPHI]",default="GEPHI")
parser.add_argument("--fastUF",help="Fast Unfolding folder [FAST_UNFOLD]",default="FAST_UNFOLD")

#  Flags
parser.add_argument("--noFastUF",help="Don't try to run Fast Unfolding Clustering, but complete everything else [False] -- Flag",action="store_true")
parser.add_argument("--usePearson",help="Use Pearson Correlation instead of Spearman [False] -- Flag",action="store_true")
parser.add_argument("--keepBigFA",help="Keep ALL big HDF5 files after finishing [False] -- Flag",action="store_true")
parser.add_argument("--keepBigF",help="Keep median correlations HDF5 files after finishing [False] -- Flag",action="store_true")
parser.add_argument("--ignoreDuplicates",help="Ignore correlation duplicates when cutting top --edgePG genes [False] -- Flag, faster if set",action="store_true")
parser.add_argument("--singleCorr",help="Output individual gene correlation files -- Warning this generates 1 file per gene in final correlation matrix. [False] -- Flag",action="store_true")
parser.add_argument("--singleCorrL",help="Output individual gene correlation files -- limited to those in --singleCorrListF [False] -- Flag",action="store_true")

#  Single corr related
parser.add_argument("--singleCorrListF",help="If --singleCorrL is set then create single correlation files for all those in this file (must be within --workFolder)",default="corrGenes.txt")

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
args.metaInf = os.path.join(args.dataF,args.metaInf)
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
	"Meta info file describing expression files [-m, --metaInf] = " + args.metaInf,
	"Separator used in expression file [-s, --fileSep] = " + args.fileSep,
	"Retain gene fraction [-f, --retainF] = " + str(args.retainF),
	"Fraction of expression files gene required in to be retained [-g, --geneFrac] = " + str(args.geneFrac),
	"Edges to retain per gene [-e, --edgePG] = " + str(args.edgePG),
	"\n#  Not changed often"
	"Root output folder [--outF] = " + args.outF,
	"Correlation matrix folder [--corrMatF] = " + args.corrMatF,
	"Single Gene Correlation files folder [--corrMatFS] = " + args.corrMatFS,
	"Gephi files folder [--gephiF] = " + args.gephiF,
	"Fast unfolding folder [--fastUF] = " + args.fastUF,
	"\n#  Flags",
	"Don't run Fast Unfolding clustering [--noFastUF] = " + str(args.noFastUF),
	"Use Pearson Correlation [--usePearson] = " + str(args.usePearson),
	"Keep big HDF5 files after run [--keepBigFA] = " + str(args.keepBigFA),
	"Keep median correlations HDF5 files after run [--keepBigF] = " + str(args.keepBigF),
	"Ignore correlation duplicates when cutting top --edgePG genes [--ignoreDuplicates] = " + str(args.ignoreDuplicates),
	"Output individual gene correlation files [--singleCorr] = " + str(args.singleCorr),
	"Output individual gene correlation files for select list (--singleCorrListF) [--singleCorrL] = " + str(args.singleCorrL),
	"\n#  Single gene correlation options",
	"List of genes to process if --singleCorrL [--singleCorrListF]:\t" + str(args.singleCorrListF),
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


def makeSafeFilename(inputFilename):
	#  Declare safe characters here
	safechars = string.letters + string.digits + "~-_.@#"
	try:
		return filter(lambda c: c in safechars, inputFilename)
	except:
		return ""
	pass


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


def listToPercentiles(x,inverse=False):
	data = stats.rankdata(x, "average")/float(len(x))

	if inverse:
		return (1-data) * 100
	else:
		return data * 100


def dataSpread(x,fixNeg=True):
	"""
	Returns the min, Q1 (25%), median (Q2), Q3 (75%), max, IQR, Quartile Coefficient of Dispersion and IQR/Median (CV like)

	As QCOD can't cope with negative numbers if fixNeg==True will solve the problem by simply adding the abs(min) value to
	all values.  So new minimum will be zero.
	"""

	if fixNeg:
		if min(x) < 0:
			x = np.array(x)
			x = x + abs(min(x))
			
	q1 = float(np.percentile(x,25,interpolation="lower"))
	q2 = np.percentile(x,50)
	q3 = float(np.percentile(x,75,interpolation="higher"))

	if (q2 == 0) or ((q3+q1) == 0):
		return min(x), q1, q2, q3, max(x), (q3-q1), 0,0
	else:
		return min(x), q1, q2, q3, max(x), (q3-q1), (q3-q1)/(q3+q1), (q3-q1)/q2


def mad(arr):
	""" Median Absolute Deviation: a "Robust" version of standard deviation.
		Indices variabililty of the sample.
	"""
	arr = np.array(arr)
	med = np.median(arr)
	return np.median(np.abs(arr - med))


def loadCorrGeneList(workF,listFN,defaultG="IRF4"):

	print("\nLoading list of genes for single correlation file output")
	fileP = os.path.join(workF,listFN)

	if not os.path.exists(fileP):
		print("\t\t",fileP,"does not exist, will create file and add default gene (",defaultG,")")
		outF = open(fileP,"w")
		print(defaultG,file=outF)
		outF.close()

		corrG = {defaultG:1}

	else:
		corrG = {}
		for line in open(fileP):
			line = line.rstrip()

			corrG[line] = 1

	print("\tTotal genes:",len(corrG))

	return corrG


def generateGeneMeta(outF,fileN,geneMetaInfo,genesNP,npA):
	"""
	Generate meta information for gene
	"""
	###########################################################################################
	def expToPercentiles(expToPerc,expA):
		return [expToPerc[e] for e in expA]

	###########################################################################################
	print("\t\t# Generate gene meta information : ",end="")
	
	expressionVals = npA.flatten()

	#  Convert from expression values to percentiles
	expPercentiles = listToPercentiles(expressionVals,inverse=False)

	#  Calculate mapping from expression to percentile
	expToPerc = {}
	for i,e in enumerate(expressionVals):
		expToPerc[e] = expPercentiles[i]
	
	outF = open(os.path.join(outF,fileN + "_metaInf.txt"),"w")
	for i,gene in enumerate(genesNP):
		#  Convert from expression values to percentiles, so can work out MAD of percentiles
		genePercentiles = expToPercentiles(expToPerc,npA[i])
		
		#  min, Q1 (25%), median (Q2), Q3 (75%), max, IQR, Quartile Coefficient of Dispersion and IQR/Median (CV like)
		try:
			minE,q1E,q2E,q3E,maxE,iqrE,qcodE,iqrME = dataSpread(npA[i])
		except:
			print("Issue with :",gene,npA[i])
			sys.exit()

		medianPercentiles = np.median(genePercentiles)
		print(gene,q2E,qcodE,medianPercentiles,sep="\t",file=outF)

		geneMetaInfo[gene].append([q2E,qcodE,medianPercentiles])

	outF.close()
	print("done")


def mergeMetaInfo(geneMetaInfo):
	print("\n\tMerging gene meta information")
	medianMetaInfo = {}

	for gene in geneMetaInfo:
		medians, qcods, percentiles = zip(*geneMetaInfo[gene])
		#  Store Median_QCOD(varWithin), MedianPercentile, MADpercentile(VarAcross)
		medianMetaInfo[gene] = [np.median(qcods),np.median(percentiles),mad(percentiles)]

	return medianMetaInfo


def loadMeta(metaFile,dataF,splitBy="\t",headerL=1):

	print("\n\nLoad meta-file (",os.path.basename(metaFile),")",sep="")
	if not os.path.exists(metaFile):
		print("\t\t# Meta file (",metaFile,") does not exist!",sep="")
		sys.exit()

	header = headerL
	fileInfo = {}
	for line in open(metaFile):
		cols = line.rstrip().split(splitBy)

		if header:
			header -= 1
			continue
		
		if line.rstrip() == "":
			continue

		totalPath = os.path.join(dataF,cols[0])
		if not os.path.exists(totalPath):
			print("\t\t# File (",totalPath,") does not exist!, won't add to fileInfo",sep="")
		else:
			try:
				fileInfo[totalPath] = int(cols[1])
			except:
				print("Meta file line (",line.rstrip(),") is not formed properly, skipping")

	print("\tLoaded information on:",len(fileInfo),"files")
	return fileInfo


def getGenes(metaInf,splitBy="\t",geneFrac=1/3.0,retainF=0.8):
	"""
	First pass over expression data-files to find the set of genes we will be
	working with after filtering each data-set by retainF and then only keeping
	genes that are present in geneFrac
	"""
	natSort = make_key_naturalSort()

	geneCount = defaultdict(int)

	print("\nLoad information on genes present per expression array:")
	for fileP in sorted(metaInf,key=natSort):
		fileN = os.path.basename(fileP)

		print("\t",fileN,end="")

		headerL = metaInf[fileP]
		
		#############################
		# Read in data file
		seenGenes = {}
		genes = []
		expression = []

		print("\tReading expression data file:",fileN)

		for i, line in enumerate(open(fileP)):
			cols = line.rstrip().split(splitBy)

			if headerL:
				headerL -= 1
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
		print("\t\t# Generate Numpy matrices")
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

		#  Keep track of seen genes
		for gene in genesNP:
			geneCount[gene] += 1

	print("\n\tTotal unique genes/identifiers:",len(geneCount))

	requiredFileNum = int(round(geneFrac * len(metaInf),0))
	print("\tRetaining genes present in ",round(geneFrac,2)," (>=",requiredFileNum,") of files : ",end="",sep="")

	genesToKeep = {}
	for g in geneCount:
		if geneCount[g] >= requiredFileNum:
			genesToKeep[g] = 1

	print(len(genesToKeep))

	return genesToKeep,requiredFileNum


def loadHDF5files(outF,hdf5Paths):

	hdf5Files = {}
	for fileN in hdf5Paths:
		hdf5Files[fileN] = h5py.File(fileN,'r')

	return hdf5Files


def generateMasks(sortedKeepGenes,genesPerHDF5,hdf5paths):

	print("\tGenerating index mappings and masks to speed up combining correlations")
	
	masks = np.zeros((len(hdf5paths),len(sortedKeepGenes)))
	
	for i,g in enumerate(sortedKeepGenes):
		for j,hd5 in enumerate(hdf5paths):
			if g not in genesPerHDF5[hd5]:
				masks[j][i] = 1  # Mask missing values
	return masks


def generateCorrMatrix(workF,metaInf,genesToKeep,requiredFileNum,singleCorrList,geneFrac=1/3.0,retainF=0.8,corrMatF="CORR_MATRIX",corrMatFS="CORR_MATRIX_SG",usePearson=False,corrChunk=5000,splitBy="\t",decimalP=3,printEveryDiv=10,keepBigFA=False,singleCorr=False,
	singleCorrL=False):
	"""
	Memory efficient method for generating all pairwise correlations of genes (rows) across a set of samples (columns).  Uses HDF5 files to greatly
	reduce memory useage, keeping most data residing on the hard disk.
	"""
	print("\n\nGenerate correlations for expression data:")
	natSort = make_key_naturalSort()
	sortedKeepGenes = sorted(genesToKeep,key=natSort)
	
	#  Find mapping of keep genes
	keepGeneMapping = {}
	for i, g in enumerate(sortedKeepGenes):
		keepGeneMapping[g] = i

	# Create subfolder
	outF = os.path.join(workF,corrMatF + "_GMF" + str(requiredFileNum))
	if not os.path.exists(outF):
		os.makedirs(outF)

	if singleCorr or singleCorrL:
		outFS = os.path.join(workF,corrMatFS + "_GMF" + str(requiredFileNum))
		if not os.path.exists(outF):
			os.makedirs(outF)
	
	geneMetaInfo = defaultdict(list)  # Per data-set store gene meta information
	genePosPerHDF5 = defaultdict(dict)  # Mapping of location of gene in HDF5 file
	perHDF5index = defaultdict(list)  # Mapping of gene to where they should be in final matrix
	genesPerHDF5 = defaultdict(dict)  # Keep tally of which genes are in which HDF5 file
	hdf5paths = []

	print("\n\tCalculating correlations for data file:")
	for fileP in sorted(metaInf,key=natSort):
		fileN = os.path.basename(fileP)

		print("\t",fileN,end="")

		headerL = metaInf[fileP]
		#############################
		# Read in data file
		seenGenes = {}
		genes = []
		expression = []

		print("\tReading expression data file:",fileN)

		for i, line in enumerate(open(fileP)):
			cols = line.rstrip().split(splitBy)

			if headerL:
				headerL -= 1
				continue

			#  Only retain required genes
			if cols[0] not in genesToKeep:
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
		print("\t\t# Generate Numpy matrices")
		#  Move genes to numpy array
		genesNP = np.empty((len(genes)),dtype=str)
		genesNP = np.array(genes)

		#  Store position of gene in matrix for this file and create mapping index for HDF5 files
		outMatrixHDF5 = os.path.join(outF,fileN + "_RetainF" + str(retainF) + ".h5")
		hdf5paths.append(outMatrixHDF5)

		for i,g in enumerate(genesNP):
			perHDF5index[outMatrixHDF5].append(keepGeneMapping[g])
			genePosPerHDF5[g][outMatrixHDF5] = i
			genesPerHDF5[outMatrixHDF5][g] = 1
			
		#  Create numpy array
		nrow = len(expression[0])
		ncol = len(expression)
		npA = np.zeros((ncol,nrow),dtype=np.float64)
		for i,eSet in enumerate(expression):
			npA[i] = eSet

		print("\t\t\tMarix shape:",npA.shape)

		#  Output genes
		outMatrixN = os.path.join(outF,fileN + "_RetainF" + str(retainF) + "_Genes.txt")
		np.savetxt(outMatrixN,genesNP,fmt="%s",delimiter="\t")
		
		#  Generate gene meta information
		generateGeneMeta(outF,fileN + "_RetainF" + str(retainF),geneMetaInfo,genesNP,npA)

		#######################################
		# Calculate correlations
		print("\t\t# Calculating correlations using HDF5 file to save memory")
		rowN, colN = npA.shape
		
		if os.path.exists(outMatrixHDF5):
			print("\t\tAlready exists -- skipping")
			continue
		with h5py.File(outMatrixHDF5, "w") as f:
			h5 = f.create_dataset("corr", (rowN,rowN), dtype="f8")

		# # Load into memory
		print("\t\t\tCreating HDF5 file")
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
		print("\t\t\tCalculating correlations for Chunk:")
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

		#  Remove matrices to save memory
		del npA
		del h5

	#  Calculate median gene meta information
	medianMetaInfo = mergeMetaInfo(geneMetaInfo)


	###########################################################################################
	##----------------------------Calculate Median Correlations------------------------------##
	###########################################################################################
	#  Calculate median correlation matrix
	print("\nCalculating median correlations")
	if len(hdf5paths) == 1:
		#  If we're only analysing a single data-set
		if not (singleCorr or singleCorrL):
			print("\t\tOnly single data-set analysed, skipping generating median correlations")
			#  No need to calculate median correlations, just return path to HDF5 file and genes matrix
			return outMatrixHDF5, genesNP
		else:
			print("\t\tOnly single data-set analysed, but --singeCorr/--singleCorrL so will proceed with output")

	printEvery = int(round(len(genesToKeep)/float(printEveryDiv),0))
	printE = printEvery
	tell = printE
	count = 0

	#  Create HDF5 median correlation matrix
	outMatrixHDF5 = os.path.join(outF,"#Median_RetainF" + str(retainF) + ".h5")
	with h5py.File(outMatrixHDF5, "w") as f:
		h5 = f.create_dataset("corr", (len(genesToKeep),len(genesToKeep)), dtype="f8")
	h5 = h5py.File(outMatrixHDF5,'r+')
	
	# Load HDF5 correlation files
	hdf5Files = loadHDF5files(outF,hdf5paths)
	
	#  Output genes
	outMatrixN = os.path.join(outF,"#Median_RetainF" + str(retainF) + "_Genes.txt")
	np.savetxt(outMatrixN,sortedKeepGenes,fmt="%s",delimiter="\t")

	#  Get masks
	maskMappings = generateMasks(sortedKeepGenes,genesPerHDF5,hdf5paths)
	
	print("\tCalculating median correlations (report every 1/",printEveryDiv,"th of total):",sep="")
	if singleCorr or singleCorrL:
		print("\t\tOutputting single gene correlation files:")

	for genePos, gene in enumerate(sortedKeepGenes):
		# print("\t\t",gene1)
		#  Inform user of position
		count += 1
		if count == printE:
			printE = printE + tell
			if singleCorr:
				print("\n\t\tProcessed:",count,end="\n\n")
			else:
				print("\t\t",count)
				
		rowsPerHDF5 = {}
		maskPos = []
		dataSetNames = []
		#  Grab row for gene across files
		for i,hdf5 in enumerate(hdf5paths):
			try:
				rowsPerHDF5[hdf5] = hdf5Files[hdf5]["corr"][genePosPerHDF5[gene][hdf5]]
				maskPos.append(i)
				dataSetNames.append(os.path.basename(hdf5)[:-3])
			except:
				pass

		#  Second pass
		#  Convert to numpy array
		npA = np.full((len(rowsPerHDF5),len(sortedKeepGenes)),-10,dtype=np.float64)  # Missing set to -10
		for i,hdf5 in enumerate(sorted(rowsPerHDF5,key=natSort)):
			npA[i][perHDF5index[hdf5]] = rowsPerHDF5[hdf5]  # Use indexes to place in correct location

		#  Get appropriate masks
		tempMask = []
		for i in maskPos:
			tempMask.append(maskMappings[i])

		npAMasked = ma.masked_array(npA,mask=tempMask)

		#  Generate medians
		medianRowCorr = np.copy(ma.median(npAMasked,axis=0))
		h5["corr"][genePos] = medianRowCorr

		###########################################################################################
		##-------------------------------SINGLE GENE CORR----------------------------------------##
		###########################################################################################
		def mergeCorrData(gene,corrInf,medianMetaInfo,medianCorr,missingVal=-10):
			finalCorrs = []
			dataSetCount = 0

			for corr in corrInf:
				if (corr <= missingVal).all():
					
					finalCorrs.append("")
				else:
					finalCorrs.append(str(round(corr,decimalP)))
					dataSetCount += 1

			roundedMeta = map(str,[round(origV,decimalP) for origV in medianMetaInfo[gene]])
			scaledCorr = dataSetCount ** medianCorr
			finalInfo = gene + "\t" + str(round(scaledCorr,decimalP)) + "\t" + str(dataSetCount) + "\t" + str(round(medianCorr,decimalP)) + "\t" + "\t".join(roundedMeta) + "\t" + "\t".join(finalCorrs)
			
			return scaledCorr,finalInfo
		###################################
		if singleCorr or singleCorrL:
			
			if singleCorrL:  # Only output genes in list
				if gene not in singleCorrList:
					continue

			subFolder = os.path.join(outFS,gene[0])
			singleCorrFName = os.path.join(subFolder,makeSafeFilename(gene) + "_corr_RetainF" + str(retainF) + ".txt.gz")
			if os.path.exists(singleCorrFName):
				continue

			print("\t\t\tSingle Corr:",makeSafeFilename(gene))

			#  Make subfolder
			if not os.path.exists(subFolder):
				os.makedirs(subFolder)
	
			singleCorrF = gzip.open(singleCorrFName,"wb",compresslevel=9)
			dataSetsH = "\t".join(dataSetNames)

			if usePearson:
				print("Gene\tNumDataSets^MedianPCC\tNumDataSets\tMedianPCC\tMedian_QCODexpression(VarWithin)\tMedianPercentile\tMADpercentile(VarAcross)",dataSetsH,sep="\t",file=singleCorrF)
			else:
				print("Gene\tNumDataSets^MedianRho\tNumDataSets\tMedianRho\tMedian_QCODexpression(VarWithin)\tMedianPercentile\tMADpercentile(VarAcross)",dataSetsH,sep="\t",file=singleCorrF)
			
			rankedByCorr = defaultdict(list)
			for i,g in enumerate(sortedKeepGenes):
				scaledCorr, info = mergeCorrData(g,npA.T[i],medianMetaInfo,medianRowCorr[keepGeneMapping[g]])
				rankedByCorr[scaledCorr].append(info)

			#  Rank by scaledCorr
			for sCorr in sorted(rankedByCorr,reverse=True):
				for info in rankedByCorr[sCorr]:
					print(info,file=singleCorrF)
			singleCorrF.close()

		###########################################################################################

	# # Write output out for debugging
	# finalCorr = np.copy(h5["corr"][:])
	# outMatrixN = os.path.join(outF,"#Median_RetainF" + str(retainF) + "_Corr")
	# stringFormat = "%." + str(decimalP) + "f"
	# np.savetxt(outMatrixN + ".txt",finalCorr,fmt=stringFormat,delimiter="\t")
	###########################################################################################
	###########################################################################################

	#  Remove all single HDF5 files unless requested to keep them
	if not keepBigFA:
		del hdf5Files
		for hdf5P in hdf5paths:
			os.remove(hdf5P)

	# #  Return path to HDF5 file and genes matrix
	return outMatrixHDF5, sortedKeepGenes


def reduceEdges(workF,dataF,gephiF,corrh5,genesM,retainF=0.8,edgePG=3,printEveryDiv=10,corrMatF=None,keepBigFA=False,keepBigF=False,ignoreDuplicates=False):
	"""
	Reduce edges in correlation matrix, only retaining edgePG maximum correlated genes per row
	"""
	
	def bottomToZero(npA,n=1):
		"""
		Set everything below n to zero
		"""
		topI = np.argpartition(npA,-n)
		npA[topI[:-n]] = 0
		return npA

	def bottomToZeroWithDuplicates(npA,n=1):
		"""
		Set everything below n to zero,
		but deal with duplicates
		"""
		unique = np.unique(npA)
		topIunique = np.argpartition(unique,-n)[-n:]

		toKeep = []
		for val in unique[topIunique]:
			toKeep.extend(np.where(npA == val)[0])

		#  Mask and reverse
		mask = np.ones(len(npA),np.bool)
		mask[toKeep] = 0
		npA[mask] = 0

		return npA

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

		if ignoreDuplicates:
			h5["corr"][i] = bottomToZero(row,edgePG + 1)
		else:
			h5["corr"][i] = bottomToZeroWithDuplicates(row,edgePG + 1)

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

	if not (keepBigFA or keepBigF):
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

	if args.singleCorrL:
		singleCorrList = loadCorrGeneList(args.workFolder,args.singleCorrListF)
	else:
		singleCorrList = None

	#  Load meta-file detailing list of files to process
	metaInf = loadMeta(args.metaInf,args.dataF)

	#  Find genes present in each data-set
	genesToKeep, requireFileNum = getGenes(metaInf,splitBy=args.fileSep,geneFrac=args.geneFrac,retainF=args.retainF)

	#  Generate correlation HDF5 files
	corrh5,genesM = generateCorrMatrix(OUT_FOLDER,metaInf,genesToKeep,requireFileNum,singleCorrList,geneFrac=args.geneFrac,retainF=args.retainF,corrMatF=args.corrMatF,corrMatFS=args.corrMatFS,usePearson=args.usePearson,corrChunk=args.corrChunk,
		splitBy=args.fileSep,keepBigFA=args.keepBigFA,singleCorr=args.singleCorr,singleCorrL=args.singleCorrL)

	#  Reduce edges
	edgesFP,nodesFP = reduceEdges(OUT_FOLDER,args.dataF,args.gephiF,corrh5,genesM,retainF=args.retainF,edgePG=args.edgePG,corrMatF=args.corrMatF + "_GMF" + str(requireFileNum),keepBigFA=args.keepBigFA,keepBigF=args.keepBigF,ignoreDuplicates=args.ignoreDuplicates)

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

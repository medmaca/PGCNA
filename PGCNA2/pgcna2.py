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

#  PGCNA2 changes  #
This is a modification of the original pgcna-multi.py that uses the
improved Leidenalg community detection method (https://arxiv.org/abs/1810.08473)
over the Louvain community detection method used originally.

"""
__author__="Matthew Care"
__version__ = "2.0.0"

import sys
import os
import argparse
import shutil
import string
from datetime import datetime
from collections import defaultdict
import re
import gzip
import random
import numpy as np
import numpy.ma as ma
import scipy.stats as stats
import h5py

##----------------------------------------Create argparse argments----------------------------------------##
parser = argparse.ArgumentParser()
#  Required
parser.add_argument("-w", "--workFolder", help="Work Folder [REQUIRED]", default=None)
parser.add_argument("-d", "--dataF", help="Expression data folder path [REQUIRED]", default=None)

#  Optional
parser.add_argument("-m", "--metaInf", help="File containing information about data files (must be in same folder as --dataF) [#FileInfo.txt]", default="#FileInfo.txt")
parser.add_argument("-s", "--fileSep", help="Separator used in expression file(s) [\t]", default="\t")
parser.add_argument("-f", "--retainF", help="Retain gene fraction -- keeping most variant genes [0.8]", type=float, default=0.8)
parser.add_argument("-g", "--geneFrac", help="Fraction of files a gene needs to be present in to be included in median corr matrix [1/3]", type=float, default=1/3.0)
parser.add_argument("-e", "--edgePG", help="Edges to keep per gene [3] -- Highly recommend leaving as default", type=int, default=3)
parser.add_argument("-r", "--roundL", help="Decimal places to round to before edge reduction [3]", type=int, default=3)

#  Not often changed
parser.add_argument("--outF", help="Root output folder [PGCNA]", default="PGCNA")
parser.add_argument("--corrMatF", help="Correlation matrix folder [CORR_MATRIX]", default="CORR_MATRIX")
parser.add_argument("--corrMatFS", help="Folder for single gene correlation files [CORR_MATRIX_SG]", default="CORR_MATRIX_SG")
parser.add_argument("--gephiF", help="Folder to store files for Gephi [GEPHI]", default="GEPHI")

#  Flags
parser.add_argument("--noLeidenalg", help="Don't try to run Leidenalg clustering, but complete everything else [False] -- Flag", action="store_true")
parser.add_argument("--usePearson", help="Use Pearson Correlation instead of Spearman [False] -- Flag", action="store_true")
parser.add_argument("--keepBigFA", help="Keep ALL big HDF5 files after finishing [False] -- Flag", action="store_true")
parser.add_argument("--keepBigF", help="Keep median correlations HDF5 files after finishing [False] -- Flag", action="store_true")
parser.add_argument("--ignoreDuplicates", help="Ignore correlation duplicates when cutting top --edgePG genes [False] -- Flag, faster if set", action="store_true")
parser.add_argument("--singleCorr", help="Output individual gene correlation files -- Warning this generates 1 file per gene in final correlation matrix. [False] -- Flag", action="store_true")
parser.add_argument("--singleCorrL", help="Output individual gene correlation files -- limited to those in --singleCorrListF [False] -- Flag", action="store_true")

#  Single corr related
parser.add_argument("--singleCorrListF", help="If --singleCorrL is set then create single correlation files for all those in this file (must be within --workFolder)", default="corrGenes.txt")

#  Correlation related
parser.add_argument("--corrChunk", dest="corrChunk", help="Size of chunk (rows) to split correlation problem over [5000] -- Higher will speed up correlation calculation at cost of RAM", default=5000, type=float)

#  Leidenalg (community detection) specific
parser.add_argument("-n", "--laNumber", dest="laRunNum", help="Number of times to run Leidenalg [100]", default=100, type=int)
parser.add_argument("-b", dest="laBestPerc", help='Copy top [10]  %% of clusterings into lBaseFold/BEST folder', default=10, type=int)
parser.add_argument("--lBaseFold", dest="lBaseFold", help="Leidenalg base folder [LEIDENALG]", default="LEIDENALG")
parser.add_argument("--lClustTxtFold", dest="lClustTxtFold", help="Leidenalg Clusters text folders [ClustersTxt]", default="ClustersTxt")
parser.add_argument("--lClustListFold", dest="lClustListFold", help="Leidenalg Clusters module list folder [ClustersLists]", default="ClustersLists")


###########################################################################################
args = parser.parse_args()
if not args.workFolder or not args.dataF:
	print("\n\nNeed to specifiy REQUIRED variables see help (-h)")
	sys.exit()
###########################################################################################
OUT_FOLDER = os.path.join(args.workFolder, args.outF, "EPG" + str(args.edgePG))
if not os.path.exists(OUT_FOLDER):
	os.makedirs(OUT_FOLDER)

###########################################################################################
args.metaInf = os.path.join(args.dataF, args.metaInf)
###########################################################################################
if not args.noLeidenalg:
	try:
		import igraph as ig
	except ModuleNotFoundError:
		print("\nCould not import igraph module need to make sure it is installed (e.g. conda install -c vtraag python-igraph ).  To skip clustering step use --noLeidenalg flag, see readme for more information\n")
		sys.exit()

	try:
		import leidenalg as la
	except ModuleNotFoundError:
		print("\nCould not import leidenalg module need to make sure it is installed (e.g. conda install -c vtraag leidenalg ).  To skip clustering step use --noLeidenalg flag, see readme for more information\n")
		sys.exit()


###########################################################################################
#  Tidy input
#  Decode fileSep for passed "\t"
args.fileSep = args.fileSep.encode().decode('unicode_escape')  # string-escape in py2
args.corrChunk = int(round(args.corrChunk, 0))
args.laRunNum = int(round(args.laRunNum, 0))

##----------------------------------------LOGGER----------------------------------------##


class multicaster(object):
	def __init__(self, filelist):
		self.filelist = filelist

	def write(self, str):
		for f in self.filelist:
			f.write(str)


def concatenateAsString(joinWith, *args):
	temp = [str(x) for x in args]
	return joinWith.join(temp)

print("##----------------------------------------", str(datetime.now()), "----------------------------------------##", sep="")
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
	"Decimal places to round to before edge reduction [-r, --roundL] = " + str(args.roundL),
	
	"\n#  Main Folders"
	"Root output folder [--outF] = " + args.outF,
	"Correlation matrix folder [--corrMatF] = " + args.corrMatF,
	"Single Gene Correlation files folder [--corrMatFS] = " + args.corrMatFS,
	"Gephi files folder [--gephiF] = " + args.gephiF,

	"\n#  Flags",
	"Don't run Fast Unfolding clustering [--noLeidenalg] = " + str(args.noLeidenalg),
	"Use Pearson Correlation [--usePearson] = " + str(args.usePearson),
	"Keep all big HDF5 files after run [--keepBigFA] = " + str(args.keepBigFA),
	"Keep median correlations HDF5 files after run [--keepBigF] = " + str(args.keepBigF),
	"Ignore correlation duplicates when cutting top --edgePG genes [--ignoreDuplicates] = " + str(args.ignoreDuplicates),
	"Output individual gene correlation files [--singleCorr] = " + str(args.singleCorr),
	"Output individual gene correlation files for select list (--singleCorrListF) [--singleCorrL] = " + str(args.singleCorrL),

	"\n#  Single gene correlation options",
	"List of genes to process if --singleCorrL [--singleCorrListF]:\t" + str(args.singleCorrListF),

	"\n#  Correlation Options",
	"Chunk size (rows) to split correlation over [--corrChunk]:\t" + str(args.corrChunk),

	"\n#  Leidenalg Specific",
	"Run number [-n, --laNumber]:\t" + str(args.laRunNum),
	"Copy top % of clusterings into *_BEST folder [-b, --laBestPerc]:\t" + str(args.laBestPerc),
	"Base folder [--lBaseFold] = " + str(args.lBaseFold),
	"Clusters text folder [--lClustTxtFold] = " + str(args.lClustTxtFold),
	"Clusters List folder [--lClustListFold] = " + str(args.lClustListFold),

)

settingsF = open(os.path.join(OUT_FOLDER, str(datetime.now()).replace(" ", "-").replace(":", ".")[:-7] + "_CorrelationSettingsInfo.txt"), "w")
print(settingInf, file=settingsF)
settingsF.close()
###########################################################################################


##----------------------------------------Methods----------------------------------------##


def makeSafeFilename(inputFilename):
	#  Declare safe characters here
	safechars = string.ascii_letters + string.digits + "~-_.@#"  # string.letters in py2
	try:
		return "".join(list(filter(lambda c: c in safechars, inputFilename)))
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


def listToPercentiles(x, inverse=False):
	data = stats.rankdata(x, "average")/float(len(x))

	if inverse:
		return (1-data) * 100
	else:
		return data * 100


def dataSpread(x):
	"""
	Returns the min, Q1 (25%), median (Q2), Q3 (75%), max, IQR, Quartile Coefficient of Dispersion and IQR/Median (CV like)
	"""

	q1 = float(np.percentile(x, 25, interpolation="lower"))
	q2 = np.percentile(x, 50)
	q3 = float(np.percentile(x, 75, interpolation="higher"))

	if (q2 == 0) or ((q3+q1) == 0):
		return min(x), q1, q2, q3, max(x), abs(q3-q1), 0, 0
	else:
		return min(x), q1, q2, q3, max(x), abs(q3-q1), abs((q3-q1)/(q3+q1)), abs((q3-q1)/q2)


def mad(arr):
	""" Median Absolute Deviation: a "Robust" version of standard deviation.
		Indices variabililty of the sample.
	"""
	arr = np.array(arr)
	med = np.median(arr)
	return np.median(np.abs(arr - med))


def loadCorrGeneList(workF, listFN, defaultG="IRF4"):

	print("\nLoading list of genes for single correlation file output")
	fileP = os.path.join(workF, listFN)

	if not os.path.exists(fileP):
		print("\t\t", fileP, "does not exist, will create file and add default gene (", defaultG, ")")
		outF = open(fileP, "w")
		print(defaultG, file=outF)
		outF.close()

		corrG = {defaultG: 1}

	else:
		corrG = {}
		for line in open(fileP):
			line = line.rstrip()

			corrG[line] = 1

	print("\tTotal genes:", len(corrG))

	return corrG


def generateGeneMeta(outF, fileN, geneMetaInfo, genesNP, npA):
	"""
	Generate meta information for gene
	"""
	###########################################################################################
	def expToPercentiles(expToPerc, expA):
		return [expToPerc[e] for e in expA]

	###########################################################################################
	print("\t\t# Generate gene meta information : ", end="")

	expressionVals = npA.flatten()

	#  Convert from expression values to percentiles
	expPercentiles = listToPercentiles(expressionVals, inverse=False)

	#  Calculate mapping from expression to percentile
	expToPerc = {}
	for i, e in enumerate(expressionVals):
		expToPerc[e] = expPercentiles[i]

	outF = open(os.path.join(outF, fileN + "_metaInf.txt"), "w")
	for i, gene in enumerate(genesNP):
		#  Convert from expression values to percentiles, so can work out MAD of percentiles
		genePercentiles = expToPercentiles(expToPerc, npA[i])

		#  min, Q1 (25%), median (Q2), Q3 (75%), max, IQR, Quartile Coefficient of Dispersion and IQR/Median (CV like)
		try:
			minE, q1E, q2E, q3E, maxE, iqrE, qcodE, iqrME = dataSpread(npA[i])
		except:
			print("Issue with :", gene, npA[i])
			sys.exit()

		medianPercentiles = np.median(genePercentiles)
		print(gene, q2E, qcodE, medianPercentiles, sep="\t", file=outF)

		geneMetaInfo[gene].append([q2E, qcodE, medianPercentiles])

	outF.close()
	print("done")


def mergeMetaInfo(geneMetaInfo):
	print("\n\tMerging gene meta information")
	medianMetaInfo = {}

	for gene in geneMetaInfo:
		medians, qcods, percentiles = zip(*geneMetaInfo[gene])
		#  Store Median_QCOD(varWithin), MedianPercentile, MADpercentile(VarAcross)
		medianMetaInfo[gene] = [np.median(qcods), np.median(percentiles), mad(percentiles)]

	return medianMetaInfo


def loadMeta(metaFile, dataF, splitBy="\t", headerL=1):

	print("\n\nLoad meta-file (", os.path.basename(metaFile), ")", sep="")
	if not os.path.exists(metaFile):
		print("\t\t# Meta file (", metaFile, ") does not exist!", sep="")
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

		totalPath = os.path.join(dataF, cols[0])
		if not os.path.exists(totalPath):
			print("\t\t# File (", totalPath, ") does not exist!, won't add to fileInfo", sep="")
		else:
			try:
				fileInfo[totalPath] = int(cols[1])
			except:
				print("Meta file line (", line.rstrip(), ") is not formed properly, skipping")

	print("\tLoaded information on:", len(fileInfo), "files")
	return fileInfo


def getGenes(metaInf, splitBy="\t", geneFrac=1/3.0, retainF=0.8):
	"""
	First pass over expression data-files to find the set of genes we will be
	working with after filtering each data-set by retainF and then only keeping
	genes that are present in geneFrac
	"""
	natSort = make_key_naturalSort()

	geneCount = defaultdict(int)

	print("\nLoad information on genes present per expression array:")
	for fileP in sorted(metaInf, key=natSort):
		fileN = os.path.basename(fileP)

		print("\t", fileN, end="")

		headerL = metaInf[fileP]
		#############################
		# Read in data file
		seenGenes = {}
		genes = []
		expression = []

		print("\tReading expression data file:", fileN)

		for i, line in enumerate(open(fileP)):
			cols = line.rstrip().replace('"', '').replace("'", '').split(splitBy)

			if headerL:
				headerL -= 1
				continue

			genes.append(cols[0])
			if not cols[0] in seenGenes:
				seenGenes[cols[0]] = 1
			else:
				print("\n\nERROR: Duplicate gene (", cols[0], ") in expression data (", fileN, "), expects unique identifiers!\nPlease remove duplicates and re-run\n", sep="")
				sys.exit()
			expression.append(cols[1:])

		print("\t\tTotal genes = ", len(genes))
		#############################

		#############################
		print("\t\t# Generate Numpy matrices")
		#  Move genes to numpy array, so can be sorted along with the rest
		genesNP = np.empty((len(genes)), dtype=str)
		genesNP = np.array(genes)

		#  Create numpy array
		nrow = len(expression[0])
		ncol = len(expression)
		npA = np.empty((ncol, nrow), dtype=object)
		for i, eSet in enumerate(expression):
			npA[i] = eSet

		#  Calculate SD for each gene
		npSD = np.std(npA.astype(np.float64), axis=1)

		#  Sort by SD
		sortI = np.argsort(npSD)[::-1]  # Sort ascending
		#  Sort matrix by index
		npA = npA[sortI, :]

		#  Sort genes by index
		genesNP = genesNP[sortI]

		#  Only retain top fract #
		#  Work out index to cut at
		cutPos = int(round(retainF * genesNP.shape[0], 0))
		print("\t\tRetain Fraction (", retainF, ") Cut at: ", cutPos, sep="")
		print("\t\t\tPre-cut shape:", npA.shape)

		#  Retain most variant
		npA = npA[:cutPos, ].astype(np.float64)
		genesNP = genesNP[:cutPos]
		print("\t\t\tPost-cut shape:", npA.shape)

		#  Keep track of seen genes
		for gene in genesNP:
			geneCount[gene] += 1

	print("\n\tTotal unique genes/identifiers:", len(geneCount))

	requiredFileNum = int(round(geneFrac * len(metaInf), 0))
	print("\tRetaining genes present in ", round(geneFrac, 2), " (>=", requiredFileNum, ") of files : ", end="", sep="")

	genesToKeep = {}
	for g in geneCount:
		if geneCount[g] >= requiredFileNum:
			genesToKeep[g] = 1

	print(len(genesToKeep))

	return genesToKeep, requiredFileNum


def loadHDF5files(outF, hdf5Paths):

	hdf5Files = {}
	for fileN in hdf5Paths:
		hdf5Files[fileN] = h5py.File(fileN, 'r')

	return hdf5Files


def generateMasks(sortedKeepGenes, genesPerHDF5, hdf5paths):

	print("\tGenerating index mappings and masks to speed up combining correlations")

	masks = np.zeros((len(hdf5paths), len(sortedKeepGenes)))
	for i, g in enumerate(sortedKeepGenes):
		for j, hd5 in enumerate(hdf5paths):
			if g not in genesPerHDF5[hd5]:
				masks[j][i] = 1  # Mask missing values
	return masks


def generateCorrMatrix(workF, metaInf, genesToKeep, requiredFileNum, singleCorrList, geneFrac=1/3.0, retainF=0.8, corrMatF="CORR_MATRIX", corrMatFS="CORR_MATRIX_SG", usePearson=False, corrChunk=5000, splitBy="\t", decimalP=3, printEveryDiv=10, keepBigFA=False, singleCorr=False, singleCorrL=False):
	"""
	Memory efficient method for generating all pairwise correlations of genes (rows) across a set of samples (columns).  Uses HDF5 files to greatly
	reduce memory useage, keeping most data residing on the hard disk.
	"""
	print("\n\nGenerate correlations for expression data:")
	natSort = make_key_naturalSort()
	sortedKeepGenes = sorted(genesToKeep, key=natSort)

	#  Find mapping of keep genes
	keepGeneMapping = {}
	for i, g in enumerate(sortedKeepGenes):
		keepGeneMapping[g] = i

	# Create subfolder
	outF = os.path.join(workF, corrMatF + "_GMF" + str(requiredFileNum))
	if not os.path.exists(outF):
		os.makedirs(outF)

	if singleCorr or singleCorrL:
		outFS = os.path.join(workF, corrMatFS + "_GMF" + str(requiredFileNum))
		if not os.path.exists(outF):
			os.makedirs(outF)

	geneMetaInfo = defaultdict(list)  # Per data-set store gene meta information
	genePosPerHDF5 = defaultdict(dict)  # Mapping of location of gene in HDF5 file
	perHDF5index = defaultdict(list)  # Mapping of gene to where they should be in final matrix
	genesPerHDF5 = defaultdict(dict)  # Keep tally of which genes are in which HDF5 file
	hdf5paths = []

	print("\n\tCalculating correlations for data file:")
	for fileP in sorted(metaInf, key=natSort):
		fileN = os.path.basename(fileP)

		print("\t", fileN, end="")

		headerL = metaInf[fileP]
		#############################
		# Read in data file
		seenGenes = {}
		genes = []
		expression = []

		print("\tReading expression data file:", fileN)

		for i, line in enumerate(open(fileP)):
			cols = line.rstrip().replace('"', '').replace("'", '').split(splitBy)

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
				print("\n\nERROR: Duplicate gene (", cols[0], ") in expression data (", fileN, "), expects unique identifiers!\nPlease remove duplicates and re-run\n", sep="")
				sys.exit()
			expression.append(cols[1:])

		print("\t\tTotal genes = ", len(genes))
		#############################

		#############################
		print("\t\t# Generate Numpy matrices")
		#  Move genes to numpy array
		genesNP = np.empty((len(genes)), dtype=str)
		genesNP = np.array(genes)

		#  Store position of gene in matrix for this file and create mapping index for HDF5 files
		outMatrixHDF5 = os.path.join(outF, fileN + "_RetainF" + str(retainF) + ".h5")
		hdf5paths.append(outMatrixHDF5)

		for i, g in enumerate(genesNP):
			perHDF5index[outMatrixHDF5].append(keepGeneMapping[g])
			genePosPerHDF5[g][outMatrixHDF5] = i
			genesPerHDF5[outMatrixHDF5][g] = 1

		#  Create numpy array
		nrow = len(expression[0])
		ncol = len(expression)
		npA = np.zeros((ncol, nrow), dtype=np.float64)
		for i, eSet in enumerate(expression):
			npA[i] = eSet

		print("\t\t\tMarix shape:", npA.shape)

		#  Output genes
		outMatrixN = os.path.join(outF, fileN + "_RetainF" + str(retainF) + "_Genes.txt")
		np.savetxt(outMatrixN, genesNP, fmt="%s", delimiter="\t")

		#  Generate gene meta information
		generateGeneMeta(outF, fileN + "_RetainF" + str(retainF), geneMetaInfo, genesNP, npA)

		#######################################
		# Calculate correlations
		print("\t\t# Calculating correlations using HDF5 file to save memory")
		rowN, colN = npA.shape

		if os.path.exists(outMatrixHDF5):
			print("\t\tAlready exists -- skipping")
			continue
		with h5py.File(outMatrixHDF5, "w") as f:
			h5 = f.create_dataset("corr", (rowN, rowN), dtype="f8")

		# # Load into memory
		print("\t\t\tCreating HDF5 file")
		h5 = h5py.File(outMatrixHDF5, 'r+')

		if not usePearson:
			# Note this isn't a perfect Spearman, as ties are not dealt with, but it is fast, and given we're only
			# Retaining the most correlated edges will have almost no effect on the output.
			npA = npA.argsort(axis=1).argsort(axis=1).astype(np.float64)

		# subtract means from the input data
		npA -= np.mean(npA, axis=1)[:, None]

		# normalize the data
		npA /= np.sqrt(np.sum(npA*npA, axis=1))[:, None]

		# Calculate correlations per chunk
		print("\t\t\tCalculating correlations for Chunk:")
		for r in range(0, rowN, corrChunk):
			print("\t\t\t", r)
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
	outMatrixHDF5_orig = outMatrixHDF5  # Retain incase only single data-set
	outMatrixHDF5 = os.path.join(outF, "#Median_RetainF" + str(retainF) + ".h5")
	if not os.path.exists(outMatrixHDF5):
		if len(hdf5paths) == 1:
			#  If we're only analysing a single data-set
			if not (singleCorr or singleCorrL):
				print("\t\tOnly single data-set analysed, skipping generating median correlations")
				#  No need to calculate median correlations, just return path to HDF5 file and genes matrix
				return outMatrixHDF5_orig, genesNP
			else:
				print("\t\tOnly single data-set analysed, but --singeCorr/--singleCorrL so will proceed with output")

		printEvery = int(round(len(genesToKeep)/float(printEveryDiv), 0))
		printE = printEvery
		tell = printE
		count = 0

		#  Create HDF5 median correlation matrix

		with h5py.File(outMatrixHDF5, "w") as f:
			h5 = f.create_dataset("corr", (len(genesToKeep), len(genesToKeep)),dtype="f8")
		h5 = h5py.File(outMatrixHDF5, 'r+')

		# Load HDF5 correlation files
		hdf5Files = loadHDF5files(outF, hdf5paths)

		#  Output genes
		outMatrixN = os.path.join(outF, "#Median_RetainF" + str(retainF) + "_Genes.txt")
		np.savetxt(outMatrixN, sortedKeepGenes, fmt="%s", delimiter="\t")

		#  Get masks
		maskMappings = generateMasks(sortedKeepGenes, genesPerHDF5, hdf5paths)

		print("\tCalculating median correlations (report every 1/", printEveryDiv, "th of total):", sep="")
		if singleCorr or singleCorrL:
			print("\t\tOutputting single gene correlation files:")

		for genePos, gene in enumerate(sortedKeepGenes):
			# print("\t\t",gene1)
			#  Inform user of position
			count += 1
			if count == printE:
				printE = printE + tell
				if singleCorr:
					print("\n\t\tProcessed:", count, end="\n\n")
				else:
					print("\t\t", count)

			rowsPerHDF5 = {}
			maskPos = []
			dataSetNames = []
			#  Grab row for gene across files
			for i, hdf5 in enumerate(hdf5paths):
				try:
					rowsPerHDF5[hdf5] = hdf5Files[hdf5]["corr"][genePosPerHDF5[gene][hdf5]]
					maskPos.append(i)
					dataSetNames.append(os.path.basename(hdf5)[:-3])
				except:
					pass

			#  Second pass
			#  Convert to numpy array
			npA = np.full((len(rowsPerHDF5), len(sortedKeepGenes)), -10, dtype=np.float64)  # Missing set to -10
			for i, hdf5 in enumerate(sorted(rowsPerHDF5, key=natSort)):
				npA[i][perHDF5index[hdf5]] = rowsPerHDF5[hdf5]  # Use indexes to place in correct location

			#  Get appropriate masks
			tempMask = []
			for i in maskPos:
				tempMask.append(maskMappings[i])

			npAMasked = ma.masked_array(npA, mask=tempMask)

			#  Generate medians
			medianRowCorr = np.copy(ma.median(npAMasked, axis=0))
			h5["corr"][genePos] = medianRowCorr

			###########################################################################################
			##-------------------------------SINGLE GENE CORR----------------------------------------##
			###########################################################################################
			def mergeCorrData(gene, corrInf, medianMetaInfo, medianCorr, missingVal=-10):
				finalCorrs = []
				dataSetCount = 0

				for corr in corrInf:
					if (corr <= missingVal).all():
						finalCorrs.append("")
					else:
						finalCorrs.append(str(round(corr, decimalP)))
						dataSetCount += 1

				roundedMeta = map(str, [round(origV, decimalP) for origV in medianMetaInfo[gene]])
				scaledCorr = dataSetCount ** medianCorr
				finalInfo = gene + "\t" + str(round(scaledCorr, decimalP)) + "\t" + str(dataSetCount) + "\t" + str(round(medianCorr, decimalP)) + "\t" + "\t".join(roundedMeta) + "\t" + "\t".join(finalCorrs)

				return scaledCorr, finalInfo
			###################################
			if singleCorr or singleCorrL:
				if singleCorrL:  # Only output genes in list
					if gene not in singleCorrList:
						continue

				subFolder = os.path.join(outFS, gene[0])
				singleCorrFName = os.path.join(subFolder, str(makeSafeFilename(gene)) + "_corr_RetainF" + str(retainF) + ".txt.gz")
				if os.path.exists(singleCorrFName):
					continue

				print("\t\t\tSingle Corr:", str(makeSafeFilename(gene)))

				#  Make subfolder
				if not os.path.exists(subFolder):
					os.makedirs(subFolder)

				singleCorrF = gzip.open(singleCorrFName, "wt", compresslevel=9)  # "wb" in py2
				dataSetsH = "\t".join(dataSetNames)

				if usePearson:
					print("Gene\tNumDataSets^MedianPCC\tNumDataSets\tMedianPCC\tMedian_QCODexpression(VarWithin)\tMedianPercentile\tMADpercentile(VarAcross)", dataSetsH, sep="\t", file=singleCorrF)
				else:
					print("Gene\tNumDataSets^MedianRho\tNumDataSets\tMedianRho\tMedian_QCODexpression(VarWithin)\tMedianPercentile\tMADpercentile(VarAcross)", dataSetsH, sep="\t", file=singleCorrF)

				rankedByCorr = defaultdict(list)
				for i, g in enumerate(sortedKeepGenes):
					scaledCorr, info = mergeCorrData(g, npA.T[i], medianMetaInfo, medianRowCorr[keepGeneMapping[g]])
					rankedByCorr[scaledCorr].append(info)

				#  Rank by scaledCorr
				for sCorr in sorted(rankedByCorr, reverse=True):
					for info in rankedByCorr[sCorr]:
						print(info, file=singleCorrF)
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
	else:
		print("\t\tAlready exists -- skipping")

	# #  Return path to HDF5 file and genes matrix
	return outMatrixHDF5, sortedKeepGenes


def reduceEdges(workF, dataF, gephiF, corrh5, genesM, retainF=0.8, edgePG=3, printEveryDiv=10, corrMatF=None, keepBigFA=False, keepBigF=False, ignoreDuplicates=False, roundL=3):
	"""
	Reduce edges in correlation matrix, only retaining edgePG maximum correlated genes per row
	"""

	def bottomToZero(npA, n=1):
		"""
		Set everything below n to zero
		"""
		topI = np.argpartition(npA, -n)
		npA[topI[:-n]] = 0
		return npA

	def bottomToZeroWithDuplicates(npA, n=1):
		"""
		Set everything below n to zero,
		but deal with duplicates
		"""
		unique = np.unique(npA)

		uniqueGTzero = len(unique[unique > 0])
		if n > uniqueGTzero:
			#  Deal with edgePG extending into negative correlations
			n = uniqueGTzero

		topIunique = np.argpartition(unique, -n)[-n:]
		toKeep = []
		for val in unique[topIunique]:
			toKeep.extend(np.where(npA == val)[0])

		#  Mask and reverse
		mask = np.ones(len(npA), np.bool)
		mask[toKeep] = 0
		npA[mask] = 0

		return npA

	###########################################################################################
	#  Setup #
	print("\nReduces edges to (", edgePG, ") per gene:", sep="")
	fileN = os.path.basename(dataF)
	#  Create subfolder
	outF = os.path.join(workF, gephiF)
	if not os.path.exists(outF):
		os.makedirs(outF)
	nodesFP = os.path.join(outF, fileN + "_RetainF" + str(retainF) + "_EPG" + str(edgePG) + "_Nodes.tsv.gz")
	edgesFP = os.path.join(outF, fileN + "_RetainF" + str(retainF) + "_EPG" + str(edgePG) + "_Edges.tsv.gz")
	if os.path.exists(nodesFP) and os.path.exists(edgesFP):
		print("\t\tEdge reduction files already exist, skipping")
		return edgesFP, nodesFP
	###########################################################################################

	print("\tLoad HDF5 file")
	#  Load old HDF5
	h5 = h5py.File(corrh5, 'r+')
	rowN, colN = h5["corr"].shape

	printEvery = int(round(rowN/float(printEveryDiv), 0))

	printE = printEvery
	tell = printE
	count = 0

	print("\tWorking (report every 1/", printEveryDiv, "th of total):", sep="")
	for i, row in enumerate(h5["corr"]):
		#  Inform user of position
		count += 1
		if count == printE:
			printE = printE + tell
			print("\t\t", count)

		#  Before edge reduction round correlations		
		row = np.round(row,decimals=roundL)

		if ignoreDuplicates:
			h5["corr"][i] = bottomToZero(row, edgePG + 1)
		else:
			h5["corr"][i] = bottomToZeroWithDuplicates(row, edgePG + 1)

	# Write output out for debugging
	# finalCorr = np.copy(h5EPG["corr"][:])
	# outMatrixN = os.path.join(outF,fileN + "_RetainF" + str(retainF) + "_NonSym")
	# stringFormat = "%." + str(3) + "f"
	# # stringFormat2 = "%." + str(decimalPpVal) + "f"
	# np.savetxt(outMatrixN + ".txt",finalCorr,fmt=stringFormat,delimiter="\t")

	###########################################################################################
	print("\nGenerating files for Gephi network visualization tool")


	#  First output list of genes (nodes)
	print("\tCreate node file")

	nodesFile = gzip.open(nodesFP, "wt", compresslevel=9)
	print("Id\tLabel\tCluster", file=nodesFile)
	for gene in genesM:
			print(gene, gene, "NA", file=nodesFile, sep="\t")
	nodesFile.close()

	#  Second, output list of edges
	print("\tCreate edges file (report every 1/", printEveryDiv, "th of total):", sep="")
	edgesFile = gzip.open(edgesFP, "wt", compresslevel=9)
	print("Source\tTarget\tWeight\tType\tfromAltName\ttoAltName", file=edgesFile)

	# Finally output edges
	printE = printEvery
	tell = printE
	count = 0
	seen = defaultdict()
	printedLines = 0

	for i, row in enumerate(h5["corr"]):
		for j, item in enumerate(row):
			if not (i == j):  # Don't output self edges...
				geneO = genesM[i]
				geneT = genesM[j]
				if not (geneO + "-@@@-" + geneT) in seen:
					#  Output info
					if not item == 0:
						print(geneO, geneT, item, "undirected", geneO, geneT, sep="\t", file=edgesFile)
						printedLines += 1

						#  Store the fact that we've see this and it's equivalent pair
						seen[geneO + "-@@@-" + geneT] = 1
						seen[geneT + "-@@@-" + geneO] = 1

		count += 1
		if count == printEvery:
			printEvery = printEvery + tell
			print("\t\t", count, "Genes", ", Edges:", printedLines)

	edgesFile.close()
	print("\n\t\tTotal printed (", printedLines, ") edges", sep="")

	if not (keepBigFA or keepBigF):
		h5.close()
		shutil.rmtree(os.path.join(workF, corrMatF))

	return edgesFP, nodesFP


def runleidenalgF(wF, edgesFP, baseFold="LEIDENALG", outFoldTorig="ClustersTxt", outFoldLorig="ClustersLists", runNum=10, bestPerc=10, printEveryDiv=10, removeExisting=True, allFold="ALL", bestFold="BEST"):
	"""
	Convert edges into igraph (https://igraph.org/python/) format and then carry out community detection using the python Leidenalg package (https://github.com/vtraag/leidenalg).  This is an improvement
	upon the Louvain method used in pgcna.py/pgcna-multi.py that ensures that all communities are well connected within the network.
	"""
	###########################################################################################

	def convertGephiToIgraph(gephiEdgeFile, splitBy="\t", header=1, outN="igraphG.temp"):
		"""
		Convert edge file into format for import to igraph
		"""
		baseFold, fileName = os.path.split(gephiEdgeFile)
		outFileP = os.path.join(baseFold, outN)
		outFile = open(outFileP, "w")

		for line in gzip.open(gephiEdgeFile):
			line = line.decode()
			if header:
				header -= 1
				continue

			cols = line.split(splitBy)
			print(" ".join(cols[0:3]), file=outFile)

		outFile.close()
		return outFileP

	def outputModules(outFoldT, outFoldL, partition, runNum=1):
		# Natural sort
		natSort = make_key_naturalSort()

		names = partition.graph.vs["name"]
		members = partition.membership
		namePerMod = defaultdict(dict)

		for n, m in zip(names, members):
			namePerMod[int(m)+1][n] = 1

		#  Output module lists
		outFileT = open(os.path.join(outFoldT, str(runNum)+".csv"), "w")
		print("Id,Modularity Class", file=outFileT)
		for m in sorted(namePerMod):
			outFileL = open(os.path.join(outFoldL, "M" + str(m) + ".txt"), "w")
			for n in sorted(namePerMod[m], key=natSort):
				print(n, file=outFileL)
				print(n, m, sep=",", file=outFileT)
			outFileL.close()
		outFileT.close()
		return namePerMod
	###########################################################################################

	print("\nGenerating Leidenalg clusterings (n=", runNum, ")", sep="")

	#  Create output folder
	if removeExisting:
		baseL = os.path.join(wF, baseFold)
		if os.path.exists(baseL):
			shutil.rmtree(baseL)

	outFoldT = os.path.join(wF, baseFold, allFold, outFoldTorig)
	outFoldL = os.path.join(wF, baseFold, allFold, outFoldLorig)
	if not os.path.exists(outFoldT):
		os.makedirs(outFoldT)
	if not os.path.exists(outFoldL):
		os.makedirs(outFoldL)

	tempP = convertGephiToIgraph(edgesFP)
	#  Input graph
	print("\tConverting gephi edge file --> igraph edge file")
	g = ig.Graph().Read_Ncol(tempP, directed=False)
	os.remove(tempP)

	# Get weights
	graphWeights = g.es["weight"]

	infoPerClustering = []
	modInfo = []
	modInfoScores = defaultdict(list)

	printEvery = int(round(runNum/float(printEveryDiv), 0))
	printE = printEvery
	tell = printE
	count = 0

	print("\tWorking (report every 1/", printEveryDiv, "th of total):", sep="")
	for rNum in range(1, runNum + 1):
		#  Make sub folder
		outFoldLS = os.path.join(outFoldL, "Clust" + str(rNum))
		if not os.path.exists(outFoldLS):
			os.makedirs(outFoldLS)

		seed = random.randint(0, 1e9)  # Generate random seed as leidenalg doesn't appear to do this (even though it says it does)
		partition = la.find_partition(g, la.ModularityVertexPartition, weights=graphWeights, n_iterations=-1, seed=seed)  # n_iterations=-1 to converge on best local result
		modScore = partition.modularity

		modInfo.append([modScore, partition.sizes()])
		modInfoScores[modScore].append(rNum)
		#  Output cluster membership and store for generating co-asociation matrix
		infoPerClustering.append(outputModules(outFoldT, outFoldLS, partition, runNum=rNum))
		#  Inform user of position
		count += 1
		if count == printE:
			printE = printE + tell
			#  Give some feedback about quality of results as working
			print("\t\t#", count, " ModScore:", round(modScore, 3), " ModNum:", len(partition.sizes()), " ModSizes:", partition.sizes())

	#  Print out information on mod scores/sizes
	modInfoF = open(os.path.join(wF, baseFold, "moduleInfoAll.txt"), "w")
	print("Mod#\tModularityScore\tAvgModSize\tModuleNum\tModuleSizes", file=modInfoF)
	for i,  inf in enumerate(modInfo):
		print(i+1, inf[0], sum(inf[1])/float(len(inf[1])), len(inf[1]), inf[1], sep="\t", file=modInfoF)
	modInfoF.close()

	if runNum > 1:
		###########################################################################################
		#  Copy best results (based on modularity score) to bestFold folder.
		#  Create output folder
		outFoldTB = os.path.join(wF, baseFold, bestFold, outFoldTorig)
		outFoldLB = os.path.join(wF, baseFold, bestFold, outFoldLorig)
		if not os.path.exists(outFoldTB):
			os.makedirs(outFoldTB)
		if not os.path.exists(outFoldLB):
			os.makedirs(outFoldLB)

		infoPerClusteringBest = []
		copyNumber = int(round(runNum * (bestPerc/100.0), 0))
		if copyNumber == 0:
			copyNumber = 1

		print("\tWill copy Best ", bestPerc, "% (n=", copyNumber, ") results to ", os.path.join(baseFold, bestFold), sep="")
		copyNum = 0

		modInfoF = open(os.path.join(wF, baseFold, "moduleInfoBest.txt"), "w")
		print("Mod#\tModularityScore\tAvgModSize\tModuleNum\tModuleSizes", file=modInfoF)
		for mScore in sorted(modInfoScores, reverse=True):
			if copyNum >= copyNumber:
				break
			for clustNum in modInfoScores[mScore]:
				if copyNum >= copyNumber:
					break

				#  Write best module information
				inf = modInfo[clustNum - 1]
				print(clustNum, inf[0], sum(inf[1])/float(len(inf[1])), len(inf[1]), inf[1], sep="\t", file=modInfoF)

				#  Store best info per cluster info
				infoPerClusteringBest.append(infoPerClustering[clustNum - 1])

				#  Copy best results to bestFold
				textName = str(clustNum) + ".csv"
				shutil.copy(os.path.join(outFoldT, textName), os.path.join(outFoldTB, textName))
				listName = "Clust" + str(clustNum)
				shutil.copytree(os.path.join(outFoldL, listName), os.path.join(outFoldLB, listName))

				#  Increment number copied
				copyNum += 1

		modInfoF.close()
		return partition.graph.vs["name"], infoPerClusteringBest
	else:
		return partition.graph.vs["name"], infoPerClustering



###########################################################################################


def main(finishT="Finished!"):

	if args.singleCorrL:
		singleCorrList = loadCorrGeneList(args.workFolder, args.singleCorrListF)
	else:
		singleCorrList = None

	#  Load meta-file detailing list of files to process
	metaInf = loadMeta(args.metaInf, args.dataF)

	#  Find genes present in each data-set
	genesToKeep, requireFileNum = getGenes(metaInf, splitBy=args.fileSep, geneFrac=args.geneFrac, retainF=args.retainF)

	#  Generate correlation HDF5 files
	corrh5, genesM = generateCorrMatrix(OUT_FOLDER, metaInf, genesToKeep, requireFileNum, singleCorrList, geneFrac=args.geneFrac, retainF=args.retainF, corrMatF=args.corrMatF, corrMatFS=args.corrMatFS, usePearson=args.usePearson, corrChunk=args.corrChunk, splitBy=args.fileSep, 
		keepBigFA=args.keepBigFA, singleCorr=args.singleCorr, singleCorrL=args.singleCorrL)

	#  Reduce edges
	edgesFP, nodesFP = reduceEdges(OUT_FOLDER, args.dataF, args.gephiF, corrh5, genesM, retainF=args.retainF, edgePG=args.edgePG, corrMatF=args.corrMatF + "_GMF" + str(requireFileNum), keepBigFA=args.keepBigFA, keepBigF=args.keepBigF, ignoreDuplicates=args.ignoreDuplicates,roundL=args.roundL)

	#  Run Leidenalg
	if not args.noLeidenalg:
		names, infoPerClustering = runleidenalgF(OUT_FOLDER, edgesFP, runNum=args.laRunNum, bestPerc=args.laBestPerc, baseFold=args.lBaseFold, outFoldTorig=args.lClustTxtFold, outFoldLorig=args.lClustListFold)

	else:
		print("\nSkipping Leidenalg community detection and downstream output as --noLeidenalg is set")
	print(finishT)

if __name__ == '__main__':
	main()

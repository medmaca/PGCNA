#PGCNA (Parsimonious Gene Correlation Network Analysis)
##Introduction
PGCNA is a gene correlation network analysis approach that is computationally simple yet yields stable and biologically meaningful modules and allows visualisation of very large networks, showing substructure and relationships that are normally hard to see.  The parsimonious approach, retaining the 3 most correlated edges per gene, results in a vast reduction in network complexity meaning that large networks can be clustered quickly and reduced to a small output file that can be used by downstream software.

##Citation
For more details see:
[Care, M.A., Westhead, D.R., and Tooze, R.M. (2018). **Defining common principles of gene co-expression refines molecular stratification in cancer**. BioRxiv 372557](https://www.biorxiv.org/content/early/2018/07/19/372557).

Please cite this when using PGCNA.
####Paper website
PGCNA paper website: [http://pgcna.gets-it.net](http://pgcna.gets-it.net)


----------

##Requirements
Python 2.7 and the following non-standard packages : numpy, scipy (for pgcna-multi.py) and h5py.
I recommend the Anaconda distribution ([https://www.anaconda.com/download/](https://www.anaconda.com/download/)), which comes with numpy, scipy and h5py included.  This script has been tested on Windows 10 (Python 2.7.14, numpy 1.15.4, scipy 1.1.0 and h5py 2.8.0) and Linux CentOS 6.9 (python 2.7.13, numpy 1.13.1, scipy 0.19.1 and h5py 2.7.0).

###Optional, but recommended
Pgcna.py / pgcna-multi.py can be run without generating modules (**--noFastUF**) if all you require is output for network visualisation tools (e.g. Gephi/Cytoscape).  For single data-set analysis use pgcna.py, which has less non-standard package requirements and uses less memory.

####Fast unfolding of communities in large networks (Louvain)
To carry out clustering you will require the additional package **Louvain** ([Blondel, V.D., Guillaume, J.-L., Lambiotte, R., and Lefebvre, E. (2008). **Fast unfolding of communities in large networks**. J. Stat. Mech. Theory Exp. 2008, P10008.](https://arxiv.org/abs/0803.0476)).  This can be downloaded here: [https://sourceforge.net/projects/louvain/files/louvain-generic.tar.gz/download](https://sourceforge.net/projects/louvain/files/louvain-generic.tar.gz/download).  Please make sure that you've downloaded v0.3 (see package README.txt).

#####Example installation
On Linux CentOS
```
$ wget https://sourceforge.net/projects/louvain/files/louvain-generic.tar.gz
$ tar -xvzf louvain-generic.tar.gz
$ cd louvain-generic
$ make
```
Copy **louvain**, **convert** and **hierarchy** into your bin folder.

####Gephi
For visualising the output of PGCNA we highly recommend using Gephi ([https://gephi.org/](https://gephi.org/)) which is able to layout large networks very efficiently.  It has the added bonus that it includes the **louvain/FastUnfold** method to quickly visualise the modularity of the network.  See **Examples** below for details of loading the output from pgcna.py into Gephi and analysing it.


----------


##Installation

Using Mercurial: 

```
hg clone https://mcare@bitbucket.org/mcare/pythonscripts-pgcna
```

Manual download: [https://bitbucket.org/mcare/pythonscripts-pgcna/downloads/](https://bitbucket.org/mcare/pythonscripts-pgcna/downloads/)


----------


##Overview
We have endeavoured to make the process of using PGCNA as simple as possible, there are relatively few user choices.  The default parameters will give good results and most users will only want to alter the **--corrChunk** option depending on amount of available RAM (see **Common options**)

###Note on terminology
The louvain/FastUnfold method **clusters** the network, each run will produce a **clustering** of the data and each of these will contain individual **clusters/modules** (the terms are used interchangeably).  Each time louvain/FastUnfold is run, due to a random start position, it will produce a different result and though these results are highly related they differ in exact gene/module groupings.  PGCNA used the louvain/FastUnfold modularity score to rank the different clusterings of the data (number set with -n, --fuNumber parameter) and select the best (number set with -r, --fuRetain parameter) for output.

##Usage

###Single vs multiple data-set analysis
**Pgcna.py** is for processing **single data-sets**, should you want to use **multiple data-sets** to generate a median correlation file (as used in the paper) then use **pgcna-multi.py**.  Many of the parameters are shared between the two scripts, with pgcna-multi.py gaining a few additional options (highlighted below).  It is possible to analyse single data-sets with pgcna-multi.py, and this will result in exactly the same output as pgcna.py.

###Input file

####Single data-set analysis
The input file for **pgcna.py** should be an expression file that has the format of unique identifiers (genes/probes etc.) in the first column, with the remaining columns containing the expression values across the samples.  Redundancy of identifiers is not allowed and must be dealt with before processing with PGCNA.  The default is for a tab separated file with a single line of header -- however, this can be altered with the (--fileSep and --headerL parameters respectively).

####Multiple data-set analysis
When running **pgcna-multi.py** the -d/--dataF option should point to a folder containing multiple data-sets (rather than a single data-set).  These files should be the same format as those used in a single data-set analysis (i.e. unique identifiers in the first column, remaining columns containing expression values).  The data-set folder also needs an additional file (default name #FileInfo.txt; set by -m/--metaInf) that contains a tab separated list of the file names to process (1st column) and the number of header lines (to allow for meta data) in those files (2nd column), this also is expected to have a header, e.g.:
```
Filename	HeaderLines
FileName1.txt	3
FileName2.txt	1
FileName3.txt	12 
```

###Basic run
PGCNA run via **pgcna.py** script:


On windows
```
python pgcna.py -w workFolderPath -d expressionFilePath
```

On linux
```
./pgcna.py -w workFolderPath -d expressionFilePath
```

PGCNA run via **pgcna-multi.py** script:


On windows
```
python pgcna-multi.py -w workFolderPath -d expressionFilesFolderPath
```

On linux
```
./pgcna-multi.py -w workFolderPath -d expressionFilesFolderPath
```

###Common options

####Alter RAM requirements (--corrChunk)
Within the PGCNA method the step that requires the most RAM is calculating the pairwise correlations.  For a small number of genes (<=20,000) this is easily carried out in memory, but this increases dramatically with increasing gene/probe numbers: 2x10^4 = 3GB, 4x10^4 = 12GB, 6x10^4 = 27GB, 8x10^4 = 48GB, 1x10^5 = 75GB etc.  PGCNA carries out only a small portion of processing in memory, breaking the large correlation problem into chunks.  This means that with default settings (--corrChunk 5000) PGCNA can process a 2x10^5 matrix in <1GB or RAM, rather than the 300GB that would be required using memory alone.

When processing larger matrices the user can choose to increase **--corrChunk** if they want to utilise more memory for calculating the correlation.  While this will speed up the correlation step, it should be noted that as this is only one step in the PGCNA process and thus the total time saving may be minor.

Setting PGCNA to have a correlation chunk size of 20,000

```
./pgcna.py -w workFolderPath -d expressionFilePath --corrChunk 2e4
```



####Retain all genes
With default setting PGCNA will only retain the top 80% most variable genes in the expression data for analysis.  If the input expression files has been pre-filtered to remove invariant genes you may want to process all of the data, this can be acomplished using the -f/--retainF parameters:

```
./pgcna.py -w workFolderPath -d expressionFilePath --retainF 1
```

####Run without clustering
If you don't have the **louvain** package installed and wish to generate output for downstream tools you can still run PGCNA with the following command:

```
./pgcna.py -w workFolderPath -d expressionFilePath --noFastUF
```

####Changing input file separator and header size
The default input format is a tab ("\t") separated file with a single header line.  This can be changed with the **--fileSep** and **--headerL** parameters, so for a .csv file with 3 header lines:

```
./pgcna.py -w workFolderPath -d expressionFilePath --fileSep "," --headerL 3
```
Note: for pgcna-multi.py the number of header lines is set using the meta-information file (default #FileInfo.txt) rather than --headerL.  In addition all the data-sets to be processed must have the same separator, so for a folder of csv separated files:
```
./pgcna-multi.py -w workFolderPath -d expressionFilesFolderPath --fileSep ","
```

####Changing the number of retained edges
The default is to retain the top 3 most correlated edges per gene, and in our paper we show that there is no benefit from increasing this.  However, should users wish to increase this they can using the **--edgePG** parameter.  Increasing the number of edges will result in the **louvain/fastUnfold** method generating fewer clusters/module, yet these clusters will be super sets of clusters/modules produced with fewer edges.  Note: when expression data contains groups of highly correlated genes it is possible for the default 3 edges to generate orphan modules (modules disconnected from the rest of the network), in this situation increasing the number of edges should reconnect these to the total network. 

To run PGCNA retaining 5 edges 

```
./pgcna.py -w workFolderPath -d expressionFilePath --edgePG 5
```

####Louvain/FastUnfold options
By default PGCNA will run 100 clusterings of the data using louvain and will then process the best one (judged by louvain modularity score).  Users may wish to increase both the number of clusterings of the network that are carried out and the fraction retained for downstream processing.

For instance to run 10,000 clusterings and retain the top 10:

```
./pgcna.py -w workFolderPath -d expressionFilePath -n 1e4 -r 10
```

####pgcna-multi.py specific options
pgcna-multi.py merges the correlations across data-sets by calculating a gene's median correlation across all data-sets that contain that gene.  During this process if the user wants pgcna-multi.py can output a folder containing a single file per gene in the median correlation matrix.  This file shows the correlations between that single gene and all other genes across the different data-sets, allowing the user to see if the data-sets have a good level of agreement, and potentially highlighting outliers.

To output all single gene correlation files (**--singleCorr**).  Note: this **significantly** increases the time to construct the median correlation matrix:

```
./pgcna-multi.py -w workFolderPath -d expressionFilesFolderPath --singleCorr
```

If you want to limit the single gene correlations output to a subset of the total genes, you can use the **--singleCorrL** option.  This uses a file in the --workFolder (default corrGenes.txt; set by --singleCorrListF) to provide a list of genes to limit the output to.  If --singleCorrL is set and corrGenes.txt is missing it will be created and the default gene (IRF4) included.  Using --singleCorrL is highly recommended over --singleCorr if you're only interested in a small number of genes due to the significant speed increase it provides.

```
./pgcna-multi.py -w workFolderPath -d expressionFilesFolderPath --singleCorrL
```

See **Output** for information on all downstream files.

----------


##Parameters

Complete list of parameters and their default values.

Can also use built in help flag:
```
./pgcna.py -h
```

###Required
|Parameter|Description|Default Value|single/multi specific|
|---|---|---|---|
|-w, --workFolder|Base folder for all output|**Required**| |
|-d, --dataF|Expression data file path|**Required**|pgcna.py --> single file, pgcna-multi.py --> folder of files|

###Input file related
|Parameter|Description|Default Value|single/multi specific|
|---|---|---|---|
|-s, --fileSep|Separator used in expression file|"\t"| |
|--headerL|Number of header lines before expression values|1|pgcna.py only|
|-m, --metaInf|File containing information about data files to process (must be in same folder as --dataF)|#FileInfo.txt|pgcna-multi.py only|
|-f, --retainF|Retain gene fraction -- keeping most variant genes|0.8| |
|-g, --geneFrac|Fraction of expression files a gene needs to be present in to be included in median correlation matrix|1/3|pgcna-multi.py only|
|-e, --edgePG|Edges to keep per gene -- **Recommend leaving as default**|3| |

###Folders
|Parameter|Description|Default Value|single/multi specific|
|---|---|---|---|
|--outF|Root output folder|"PGCNA"| |
|--corrMatF|Correlation matrix folder -- where HDF5 files are stored|"CORR_MATRIX"| |
|--corrMatFS|Folder for single gene correlation files|"CORR_MATRIX_SG"|pgcna-multi.py only|
|--gephiF|Folder to store files for Gephi|"GEPHI"| |
|--fastUF|Folder for fast unfolding clustering output|"FAST_UNFOLD"| |

###Control usage
|Parameter|Description|Default Value|single/multi specific|
|---|---|---|---|
|--noFastUF|Flag -- don't run Fast Unfolding clustering but complete everything else|False| |
|--usePearson|Flag -- Instead of Spearman's ranked correlation coefficient calculate Pearson's|False| |
|--keepBigF|Flag -- Retain big HDF5 files after finishing, if needed for independent downstream analysis|False|For pgcna-multi.py --keepBigF retains median corr HDF5|
|--keepBigFA|Flag -- Retain **ALL** big HDF5 files after finishing, if needed for independent downstream analysis|False|pgcna-multi.py only|
|--corrChunk|Size of chunk to split correlation problem over -- higher values will speed up correlation calculation at the cost of RAM|5000| |
|--ignoreDuplicates|Flag -- Ignore correlation duplicates when cutting top --edgePG genes, faster if set but may miss genes if correlations are duplicated|False| |
|--singleCorr|Flag -- Output individual gene correlation files -- Warning this generates 1 file per gene in final correlation matrix|False|pgcna-multi.py only|
|--singleCorrL|Flag -- Output individual gene correlation files -- limited to those in --singleCorrListF|False|pgcna-multi.py only|
|--singleCorrListF|If --singleCorrL is set then create single correlation files for all those in this file (must be within --workFolder)|corrGenes.txt|pgcna-multi.py only|

###Clustering related
|Parameter|Description|Default Value|single/multi specific|
|---|---|---|---|
|-n, --fuNumber|Number of times to run **louvain/FastUnfold** method|100| |
|-r, --fuRetain|How many of the best (based on louvain modularity) clusterings to process|1| |
|--fuRenumberStart|For clusters retained (--fuRetain), what value to start numbering them from|1| |
|--tOutFold|FastUnfold Trees output folder|"Trees"| |
|--fTOutFold|FastUnfold final retained (--fuRetain) trees output folder|"TreesF"| |
|--cOutFold|FastUnfold clusters folder|"Clusters"| |
|--cOutTxtFold|FastUnfold clusters mapped back to genes|"ClustersTxt"| |
|--cOutListsFold|FastUnfold clusters mapped back to genes and split into individual text files|"ClustersLists"| |


###Output
If PGCNA is run using default settings, the root folder (-w, --workFolder) will contain the following:

* EPG3 (Edge per gene 3 folder)
	* *CorrelationSettingsInfo.txt : timestamped file containing parameters used for run.
	* CORR_MATRIX (if --keepBigF) : contains the temporary files used for calculating correlations
		* *_RetainF*.h5 : The HDF5 file of all pairwise correlations after edge reduction
		* *_RetainF*_Genes.txt : Genes that remain after -f/--retainF filtering ordered by decending standard deviation
	* CORR_MATRIX_SG (if --singleCorr/--singleCorrL) : contains gzipped single gene correlation files split into first letter subfolders
	* FAST_UNFOLD : Root folder for output from fast unfold
		* *retainF1.0_EPG3 : folder specific to an input data file
			* Clusters : Clusterings output by the FastUnfold **hierarchy** tool
			* **ClustersLists** : Contains subfolders for the -r/--fuRetain number of "best" clusterings.  Each subfolder contains genes split across the clusters/modules (M1...M*.txt)
			* **ClustersTxt** : Contains the -r/--fuRetain number of "best" clusterings as individual *.csv files
			* GeneIntMap.txt : Mapping of gene to numbers, required to get back to gene names after processing with **louvain**
			* modScores.txt : Modularity scores across all clusterings -- used to rank clusterings by and select "best" to retain (-r/--fuRetain)
			* sym_edges.bin : binary version of sym_edges.txt output by FastUnfold **convert** tool and required by **louvain**
			* sym_edges.txt : gene pairs (encoded as numbers, see GeneIntMap.txt) along with their correlation score.
			* sym_edges.weights | Edge weights, output by FastUnfold **convert** tool and require by **louvain**
			* Trees: Trees output by **louvain**
			* TreesF: Contains the -r/--fuRetain number of "best" trees.
	* **GEPHI**: contains files for processing with Gephi package ([https://gephi.org/](https://gephi.org/))
		* *_RetainF1.0_EPG3_Edges.tsv : Tab separated file of edges
		* *_RetainF1.0_EPG3_Nodes.tsv : Tab separated file of nodes

The most important folders for users are highlighted in **bold**.


##Examples

###Gephi Visualisation
To visualise the output from PGCNA in Gephi is quite straightforward (correct for Gephi V0.92):

Using data from **EPG3/GEPHI** folder

1. Open Gephi
2. Select File/"New Project"
3. Select "Data Table" (tab above screen centre)
4. Select "Import Spreadsheet"
	1. Select *_RetainF1.0_EPG3_Nodes.tsv , make sure "Import as": "Nodes table", click **Next** and then **Finish**.  Hopefully should import with no issues, select **OK**
	2. Repeat with *_RetainF1.0_EPG3_Nodes.tsv, making sure that "Import as" : "Edges table", click **Next** and then **Finish**.  Hopefully should import with no issues, select **OK**
5. Click "Graph" (tab above screen centre) -- you should now see a massive blob of nodes in the screen centre
6. Under Statistics/NetworkOverview (right of screen) select **Modularity**.  This will run the built in version of the louvain/fastUnfold method.
	1. Each run will generate a different clustering, the higher the score the better the clustering is perceived to be.
7. Select Appearance/Partition (left of screen) and select "Modularity Class" from drop-down list
	1. Click Palette button in the bottom right of this panel and select Generate. Untick "Limit number of colors" and then select **OK**.
	2. If you want you can manually change any of the resultant colours by left clicking on them and dragging the mouse around.
	3. When you're happy with the range of colours select **Apply** button.
8. To layout the network:
	1. Under **Layout** (left middle of screen) select **ForceAtlas2**
	2. For big networks set the following ("Approximate Repulsion":selected, Scaling:<1 often as small as 0.1, "LinLog mode":selected and importantly "Prevent Overlap":**not** selected).
	3. Wait until network layout has finished (You may need to alter **Scaling** if all nodes are pushed to the edge of the screen.)
	4. Once you're happy with the network set "Prevent Overlap":selected to finally prevent node overlaps.

The Gephi version of the louvain/fastUnfold is older than that used by PGCNA, so if you've run PGCNA with clustering included, you should use the results contained in the **EPG/FAST_UNFOLD/\*retainF1.0_EPG3/ClustersTxt**, assuming you've already loaded the network (see above):

1. Select "Data Table" (tab above screen centre)
2. Select "Import Spreadsheet"
	1. Within ClustersTxt choose a *.csv file to import
	2. Make sure "Import as": "Nodes table", click **Next**, make sure **Modularity Class** is set to integer then click **Finish**.
3. Follow step 7 above to set appearance with new Modularity Class data.

##Feedback and questions
If you have any queries or notice any bugs please email me at **m.a.care@leeds.ac.uk** (please include PGCNA in the subject heading).

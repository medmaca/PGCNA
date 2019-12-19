# PGCNA2 (Parsimonious Gene Correlation Network Analysis version 2)

## Introduction

PGCNA is a gene correlation network analysis approach that is computationally simple yet yields stable and biologically meaningful modules and allows visualisation of very large networks, showing substructure and relationships that are normally hard to see.  The parsimonious approach, retaining the 3 most correlated edges per gene, results in a vast reduction in network complexity meaning that large networks can be clustered quickly and reduced to a small output file that can be used by downstream software.

## Citation

For more details see:
[Care, M.A., Westhead, D.R., Tooze, R.M., 2019. Parsimonious Gene Correlation Network Analysis (PGCNA): a tool to define modular gene co-expression for refined molecular stratification in cancer. npj Syst. Biol. Appl. 5, 13. https://doi.org/10.1038/s41540-019-0090-7](https://www.nature.com/articles/s41540-019-0090-7).

Please cite this when using PGCNA.

### Paper website

PGCNA paper website: [http://pgcna.gets-it.net](http://pgcna.gets-it.net)


----------

## PGCNA2

**pgcna2.py** is a modified version of **pgcna-multi.py** that uses **Leidenalg** instead of the **Louvain** community detection algorithm.  Leidenalg has been shown to produce better clustering solutions which guarantee well-connected communities (for more details see [Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports, 9(1). https://doi.org/10.1038/s41598-019-41695-z](https://www.nature.com/articles/s41598-019-41695-z)).

## Requirements

Due to the inclusion of the leidenalg community detection method pgcna2.py requires Python 3 and the following non-standard packages : numpy, scipy, h5py (and for clustering: python-igraph and leidenalg).
We recommend the Anaconda distribution ([https://www.anaconda.com/download/](https://www.anaconda.com/download/)), which comes with numpy, scipy and h5py included and makes installing python-igraph and leidenalg simple.  This script has been tested on Windows 10 (Python 3.7.3, numpy 1.16.4, scipy 1.2.1, h5py 2.9.0, python-igraph 0.7.1 and leidenalg 0.7.0) and Linux CentOS 6.9 (Python 3.7.3, numpy 1.16.4, scipy 1.3.0, h5py 2.9.0, python-igraph 0.7.1 and leidenalg 0.7.0).

### Optional, but recommended

pgcna2.py can be run without generating modules (**--noLeidenalg**) if all you require is output for network visualisation tools (e.g. Gephi/Cytoscape).

#### Community detection with leidenalg

To carry out clustering you will require the additional python packages **python-igraph** and **leidenalg**.  

##### Example installation

On both Windows/Linux:

Presuming that you're using Anaconda python (**must be Py3 environment**) these can be installed with:
```
conda install -c vtraag python-igraph 
conda install -c vtraag leidenalg
```


#### Gephi

For visualising the output of PGCNA we highly recommend using Gephi ([https://gephi.org/](https://gephi.org/)) which is able to layout large networks very efficiently.  It has the added bonus that it includes the **louvain/FastUnfold** method to quickly visualise the modularity of the network (however, would recommend the leidenalg results over those produced by Gephi).  See **Examples** below for details of loading the output from pgcna.py into Gephi and analysing it.


----------


## Installation

Using github: 

```
git clone https://github.com/medmaca/PGCNA.git
```

Manual download: [https://github.com/medmaca/PGCNA/archive/master.zip](https://github.com/medmaca/PGCNA/archive/master.zip)


----------


## Overview

We have endeavoured to make the process of using PGCNA as simple as possible, there are relatively few user choices.  The default parameters will give good results and most users will only want to alter the **--corrChunk** option depending on amount of available RAM (see **Common options**)

### Note on terminology

The leidenalg method **clusters** the network, each run will produce a **clustering** of the data and each of these will contain individual **clusters/modules** (the terms are used interchangeably).  Each time leidenalg is run, due to a random start position, it will produce a different result and though these results are highly related they differ in exact gene/module groupings.  PGCNA used the leidenalg modularity score to rank the different clusterings of the data (number set with **-n**, **--laNumber**) and select the best (number set with **-b**, **--laBestPerc** parameter) percentage for output.

## Usage

### Data-set analysis

**pgcna2.py** can process both single and multiple data-sets (as used in the original paper). 

### Input file

When running **pgcna2** the **-d/--dataF** option should point to a folder containing one or more data-sets.  The input file(s) should be expression file(s) that have the format of unique identifiers (genes/probes etc.) in the first column, with the remaining columns containing the expression values across the samples.  Redundancy of identifiers is not allowed and must be dealt with before processing with PGCNA.  The data-set folder **also needs an additional file (default name #FileInfo.txt; set by **-m/--metaInf**)** that contains a tab separated list of the file names to process (1st column) and the number of header lines (to allow for meta data) in those files (2nd column), this also is expected to have a header, e.g.:

```
Filename	HeaderLines
FileName1.txt	3
FileName2.txt	1
FileName3.txt	12 
```

### Basic run

PGCNA run via **pgcna2.py** script:


On windows
```
python3 pgcna2.py -w workFolderPath -d expressionFilesFolderPath
```

On linux
```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath
```

### Common options

#### Alter RAM requirements (--corrChunk)

Within the PGCNA method the step that requires the most RAM is calculating the pairwise correlations.  For a small number of genes (<=20,000) this is easily carried out in memory, but this increases dramatically with increasing gene/probe numbers: 2x10^4 = 3GB, 4x10^4 = 12GB, 6x10^4 = 27GB, 8x10^4 = 48GB, 1x10^5 = 75GB etc.  PGCNA carries out only a small portion of processing in memory, breaking the large correlation problem into chunks.  This means that with default settings (--corrChunk 5000) PGCNA can process a 2x10^5 matrix in <1GB or RAM, rather than the 300GB that would be required using memory alone.

When processing larger matrices the user can choose to increase **--corrChunk** if they want to utilise more memory for calculating the correlation.  While this will speed up the correlation step, it should be noted that as this is only one step in the PGCNA process and thus the total time saving may be minor.

Setting PGCNA to have a correlation chunk size of 20,000

```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath --corrChunk 2e4
```

#### Retain all genes

With default setting PGCNA will only retain the top 80% most variable genes in the expression data for analysis.  If the input expression files has been pre-filtered to remove invariant genes you may want to process all of the data, this can be accomplished using the -f/--retainF parameters:

```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath --retainF 1
```

#### Run without clustering

If you don't have the **python-igraph/leidenalg** packages installed and wish to generate output for downstream tools you can still run PGCNA with the following command:

```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath --noLeidenalg
```

#### Changing input file separator and header size

The default input format is a tab ("\t") separated file with a single header line.  This can be changed with the **--fileSep** parameter, and the number of header lines is set using the meta-information file (default #FileInfo.txt).  In addition all the data-sets to be processed **must** have the same separator, so for a folder of csv separated files:
```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath --fileSep ","
```

#### Changing the number of retained edges

The default is to retain the top 3 most correlated edges per gene, and in our paper we show that there is no benefit from increasing this.  However, should users wish to increase this they can using the **--edgePG** parameter.  Increasing the number of edges will result in the **leidenalg** method generating fewer clusters/modules, yet these clusters will be super sets of clusters/modules produced with fewer edges.  Note: when expression data contains groups of highly correlated genes it is possible for the default 3 edges to generate orphan modules (modules disconnected from the rest of the network), in this situation increasing the number of edges should reconnect these to the total network. 

To run PGCNA retaining 5 edges 

```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath --edgePG 5
```


#### Leidenalg options

By default PGCNA will run 100 clusterings of the data using leidenalg and will retain the best 10% (judged by leidenalg modularity score).  Users may wish to increase both the number of clusterings of the network that are carried out and the fraction retained for downstream processing (copied into lBaseFold/BEST folder).

For instance to run 10,000 clusterings and retain the top 20%:

```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath -n 1e4 -b 20
```

#### Options specific to processing multiple expression files

If pgcna2.py is processing multiple expression files it merges the correlations across data-sets by calculating a gene's median correlation across all data-sets that contain that gene.  During this process if the user wants pgcna2.py can output a folder containing a single file per gene in the median correlation matrix.  This file shows the correlations between that single gene and all other genes across the different data-sets, allowing the user to see if the data-sets have a good level of agreement, and potentially highlighting outliers.

To output all single gene correlation files (**--singleCorr**).  Note: this **significantly** increases the time to construct the median correlation matrix:

```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath --singleCorr
```

If you want to limit the single gene correlations output to a subset of the total genes, you can use the **--singleCorrL** option.  This uses a file in the **--workFolder** (default corrGenes.txt; set by **--singleCorrListF**) to provide a list of genes to limit the output to.  If **--singleCorrL** is set and corrGenes.txt is missing it will be created and the default gene (IRF4) included.  Using **--singleCorrL** is highly recommended over **--singleCorr** if you're only interested in a small number of genes due to the significant speed increase it provides.

```
./pgcna2.py -w workFolderPath -d expressionFilesFolderPath --singleCorrL
```

See **Output** for information on all downstream files.

----------


## Parameters

Complete list of parameters and their default values.

Can also use built in help flag:
```
./pgcna2.py -h
```

### Required

|Parameter|Description|Default Value|
|---|---|---|
|-w, --workFolder|Base folder for all output|**Required**|
|-d, --dataF|Path to folder of expression file(s)|**Required**|

### Input file related

|Parameter|Description|Default Value|Multi specific|
|---|---|---|---|
|-s, --fileSep|Separator used in expression file|"\t"| |
|-m, --metaInf|File containing information about data files to process (must be in same folder as --dataF)|#FileInfo.txt||
|-f, --retainF|Retain gene fraction -- keeping most variant genes|0.8| |
|-g, --geneFrac|Fraction of expression files a gene needs to be present in to be included in median correlation matrix|1/3|Yes|


### Edge reduction related

|Parameter|Description|Default Value|
|---|---|---|
|-e, --edgePG|Edges to keep per gene -- **Recommend leaving as default**|3|
|-r, --roundL|Decimal places to round correlations to before edge reduction|3|


### Folders

|Parameter|Description|Default Value|Multi specific|
|---|---|---|---|
|--outF|Root output folder|"PGCNA"| |
|--corrMatF|Correlation matrix folder -- where HDF5 files are stored|"CORR_MATRIX"| |
|--corrMatFS|Folder for single gene correlation files|"CORR_MATRIX_SG"|Yes|
|--gephiF|Folder to store files for Gephi|"GEPHI"| |



### Control usage

|Parameter|Description|Default Value|Multi specific|
|---|---|---|---|
|--noLeidenalg|Flag -- don't run **leidenalg** clustering but complete everything else|False| |
|--usePearson|Flag -- Instead of Spearman's ranked correlation coefficient calculate Pearson's|False| |
|--keepBigF|Flag -- Retain big HDF5 files after finishing, if needed for independent downstream analysis|False|For multiple files --keepBigF retains median corr HDF5|
|--keepBigFA|Flag -- Retain **ALL** big HDF5 files after finishing, if needed for independent downstream analysis|False|Yes|
|--corrChunk|Size of chunk to split correlation problem over -- higher values will speed up correlation calculation at the cost of RAM|5000| |
|--ignoreDuplicates|Flag -- Ignore correlation duplicates when cutting top --edgePG genes, faster if set but may miss genes if correlations are duplicated|False| |
|--singleCorr|Flag -- Output individual gene correlation files -- Warning this generates 1 file per gene in final correlation matrix|False|Yes|
|--singleCorrL|Flag -- Output individual gene correlation files -- limited to those in --singleCorrListF|False|Yes|
|--singleCorrListF|If --singleCorrL is set then create single correlation files for all those in this file (must be within --workFolder)|corrGenes.txt|Yes|

### Clustering related

|Parameter|Description|Default Value|
|---|---|---|
|-n, --laNumber|Number of times to run **leidenalg** method|100|
|-b, --laBestPerc|Percentage of best (based on leidenalg modularity) clusterings to copy to lBaseFold/BEST |10|
|--lBaseFold|Base folder for leidenalg|"LEIDENALG"|
|--lClustTxtFold|Leidenalg Clusters text folders |ClustersTxt|
|--lClustListFold|Leidenalg Clusters module list folder |ClustersLists|



### Output

If PGCNA is run using default settings, the root folder (-w, --workFolder) will contain the following:

* **EPG3** (Edge per gene 3 folder)
	* *CorrelationSettingsInfo.txt : timestamped file containing parameters used for run.
	* **CORR_MATRIX_GMF*** (if --keepBigF) : contains the temporary files used for calculating correlations
		* *_RetainF*.h5 : The HDF5 file of all pairwise correlations after edge reduction
		* *_RetainF*_Genes.txt : Genes that remain after -f/--retainF filtering ordered by decending standard deviation
	* **CORR_MATRIX_SG_GMF*** (if --singleCorr/--singleCorrL) : contains gzipped single gene correlation files split into first letter subfolders
	* **_GEPHI_**: contains files for processing with Gephi package ([https://gephi.org/](https://gephi.org/))
		* *_RetainF1.0_EPG3_Edges.tsv.gz : Tab separated file of edges
		* *_RetainF1.0_EPG3_Nodes.tsv.gz : Tab separated file of nodes
	* **LEIDENALG** : Root folder for output from leidenalg clustering
		* **ALL** : contains all results from clustering
			* **ClustersLists** : Contains one folder per clustering each of which contains one text file per module/cluster containing that modules genes.
			* **ClustersTxt** : Contains one *.csv file per clustering each of which contains info on the gene/module pairings.
		* **BEST** : contains the best **-b/--laBestPerc** percentage results (judged by leidenalg modularity score)
			* **ClustersLists** : Contains one folder per clustering each of which contains one text file per module/cluster containing that modules genes.
			* **ClustersTxt** : Contains one *.csv file per clustering each of which contains info on the gene/module pairings.
		* moduleInfoAll.txt : Contains meta information about clusterings for the data in **ALL** folder
		* moduleInfoBest.txt : Contains meta information about clusterings for the data in **BEST** folder


## Examples

### Gephi Visualisation

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

The Gephi version of the louvain/fastUnfold is inferior to the leidenalg method used by PGCNA, so if you've run PGCNA with clustering included, you should use the results contained in the **EPG/LEIDENALG/BEST/ClustersTxt**, assuming you've already loaded the network (see above):

1. Select "Data Table" (tab above screen centre)
2. Select "Import Spreadsheet"
	1. Within ClustersTxt choose a *.csv file to import
	2. Make sure "Import as": "Nodes table", click **Next**, make sure **Modularity Class** is set to integer then click **Finish**.
3. Follow step 7 above to set appearance with new Modularity Class data.

## Troubleshooting

### Too many edges

If after running pgcna2.py the number of edges in the network is >> **-e/--edgePG** x nodes (retained genes) then there are two likely causes:

1.	The data-set has many highly related genes/probes that have identical correlations, and thus aren't reduced during edge reduction.
2.	The data-set being analysed is very sparse in nature (e.g. from single-cell RNA-seq).

In both cases setting the **--usePearson** option will potential improve the result.  This is because PGCNA uses an approximation of Spearman's correlation to reduce memory overhead and decrease run time.  However, this approximation doesn't deal with tied ranks and thus will have issues when there are blocks of highly correlated (or indeed sparse) data.

However, we would suggest that if you are suffering from this issue that you explore your data-set first to try to understand why you have so many redundant (highly correlated) genes/probes.  In addition while pgcna2.py will generate meaningful results for scRNA-seq data (with --usePearson) it has not been developed with that in mind and thus we'd recommend you seek another approach for such data.

## Feedback and questions

If you have any queries or notice any bugs please email me at **m.a.care@leeds.ac.uk** (please include PGCNA in the subject heading).

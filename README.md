TFenricher
================

About TFenricher
-------------
The transcription factor downstream annotation enricher (TFenricher) package is a bioinformatics tool to enable users to do an enrichment analysis of lists of transcription factors (TFs). 

Transcription factors (TFs) are the upstream regulators that orchestrate gene expression, and therefore a centrepiece in bioinformatics studies. A popular strategy that helps understanding the biological context of genes/proteins includes basic annotation enrichment. There are several annotation enrichment methods available, yet these methods are not well suited for analysing groups of TFs, particularly since such methods fail to include downstream processes. Here, we present TFenricher, a Python toolbox that focuses specifically at identifying gene ontologies, cellular pathways, and diseases that are overrepresented among genes that are downstream of user-defined sets of TFs. Given a set of TFs, TFenricher infers downstream genes and calculates enrichments in some of the most common databases of gene functionalities, including GO, KEGG, and Reactome. The TFenricher package enables users to search for biological context in any set of TFs and their downstream genes.

Examples of how to use TFenricher are found below, and in a jupyter notebook under ./examples/

Installation
============
The package can be used under the GNU GENERAL PUBLIC LICENSE V3

To run the default version of TFenricher, a correlation look-up table that is used in the TF-to-target mappings needs to be constructed. This is done when the setup.py file is run. We recommend that TFenricher is installed using pip:

```consol
pip install .
```

Python Version Supported/Tested
-------------------------------
- Python 3.6

Dependencies
------------
Except for the standard, built-in Python modules, the following dependencies are needed:

- [NumPy](https://www.numpy.org/)

- [pandas](https://pandas.pydata.org/)

- [Matplotlib]()

- [Scipy]()

All dependencies are included as per default in a Conda environment, 

Usage:
======
In python:
```python
# Comment
>>> import TFenricher 
```
Or from the command line:
```console
python TFenricher.py --TFs tfs.txt
```



In depth description of TFenricher
===============================
The TFenricher algorithm works in two distinct steps (Fig. 1a). First, it maps a user-defined list of TFs to putative downstream genes using lookup-tables of co-expression that comes included with the software. In detail, the expression correlation was extracted using the ARCHS4 database\cite{} and is based on data from >100k gene expression profiles, making it one of the most extensive co-expression analyses currently available. This mapping can, however, easily be replaced to a method defined by the user.

The second step takes the mapped target genes and performs enrichment analyses on gene sets annotated in, as per the choice of the user, KEGG, GO, REACTOME, the GWAS catalogue, or, alternatively, any set defined by the user. Moreover, enrichments are calculated using a Fisher’s exact test, and multiple testing correction is available using either a Bonferroni or Benjamini-Hochberg correction, or any correction provided by the user.


Contributor:
=============

 Rasmus Magnusson: Development of the package.

Current Members in the Project
------------------------------
- @rasma774

References & how to cite
======================
Add here


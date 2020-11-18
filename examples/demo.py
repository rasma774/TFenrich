#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:53:16 2020

@author: rasmus
"""

import sys
sys.path.append('../TFenrich/')
from TFenrich import TFenrich
import pandas as pd

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2020, Linköping'
__contact__ = 'rasma774@gmail.com'


TFlist = pd.read_csv('../data/gene_correlations/tfnames.txt', header=None)
ms_genes = pd.read_csv('E-GEOD-77598-analytics.tsv', sep='\t')
ms_genes = ms_genes[ms_genes.iloc[:,1].isin(TFlist[0])]
ms_genes = ms_genes.iloc[:,1][ms_genes.iloc[:,2] < 0.1]
enr = TFenrich.TFenrich(ms_genes.values, mapmethod='corr')

enr.downstream_enrich(db='REACTOME')
a = enr.pathway_enrichments.iloc[:20,:]

# now compare w just the TFs
tmp = {}
for x in ms_genes.values:
    tmp[x] = 0

tmp = pd.Series(tmp)
enr.target_genes = tmp
enr.downstream_enrich(db='REACTOME')
b = enr.pathway_enrichments.iloc[:20,:]

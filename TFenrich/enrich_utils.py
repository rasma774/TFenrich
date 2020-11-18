#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 13:26:33 2020

@author: rasmus
"""

import numpy as np
import scipy.stats as sts
import pandas as pd

from sklearn.metrics import roc_auc_score


from stat_utils import benjaminihochberg_correction

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2020, Linköping'
__contact__ = 'rasma774@gmail.com'

def _sortsets(db):
    gene_lists = {}
    f = open('../data/gene_annotations/gene_annot/c2.all.v7.1.symbols.gmt')
    for line in f:
        line = line.replace('\n', '').split('\t')
        if line[0][:len(db)] == db:
            genes = np.array(line[2:])
            gene_lists[line[0]] = genes
    f.close()
    return gene_lists

def _calc_fisher(gene_lists, genes_tmp):
      # Fisher exact test
    #               | in disease genes | not disease gene
    #----------------------------------------------------
    # in light up   |         A        |        B
    #---------------|------------------------------------
    # not light up  |         C        |        D
    #----------------------------------------------------
    #
    OR = {}
    p = {}
    for hallmark in gene_lists:
        A = np.in1d(gene_lists[hallmark], genes_tmp).sum()
        B = len(gene_lists[hallmark]) - A
        C = len(genes_tmp) - A
        D = 26408 - (A + B + C)

        OR[hallmark], p[hallmark] = sts.fisher_exact([[A, B], [C, D]], alternative='greater')
    res = pd.DataFrame([OR, p]).transpose()
    res.columns = ['OR', 'p']
    return res


def set_enrichments(gene_set, db='KEGG', FDR=0.05, ):
    """
    

    Parameters
    ----------
    gene_set : pd.array
        DESCRIPTION.
    db : str
        suggestions include, but are not limited to
        {'REACTOME', 'KEGG', 'BIOCARTA', 'PID', 'ALL'}. The default is KEGG
    FDR : float, optional
        False discovey rate acc BenjaminiHochberg. 0 < FDR < 1. The default is 0.05.
    
    Returns
    -------
    enrichment analysis .

    """
    
    gene_lists = _sortsets(db)    
    res = _calc_fisher(gene_lists, gene_set.index)
    
    index_sort = np.argsort(res.p)
    res = res.iloc[index_sort, :]
    
    res['FDR'] = benjaminihochberg_correction(res.p, FDR=FDR)
    return res

def GWAS(rankings, ismember, ngenes_thresh=100):
    """
    

    Parameters
    ----------
    rankings : TYPE
        DESCRIPTION.
    ismember : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    print('loading GWAS...')
    gwas = pd.read_csv(
        '../data/gwas/gwas.csv', 
        index_col=0,
        )
    print('Done')
    unique_dis, numb = np.unique(gwas.index, return_counts=True)
    unique_dis = unique_dis[numb >= ngenes_thresh]
    gwas = gwas.loc[unique_dis]


    rankings = rankings.copy()[rankings.index.isin(gwas.values.T[0])]    
    res ={}
    for disease in np.unique(gwas.index):
        isin = rankings.index.isin(gwas.loc[disease].values.T[0])
        p = sts.mannwhitneyu(rankings[isin].values, rankings[~isin].values)[-1]
        FE = rankings[isin].values.mean()/rankings[~isin].values.mean()
        res[disease] = [FE, p]
    res = pd.DataFrame(res).transpose()
    return res
        
def enrichr(genes):
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects.packages import importr
    #import rpy2
    from rpy2 import robjects as robj
    utils = importr('utils')
    if not rpackages.isinstalled("org.Hs.eg.db"):
        utils.install_packages("org.Hs.eg.db")
        
    if not rpackages.isinstalled("clusterProfiler"):
        utils.install_packages("clusterProfiler")
    importr("clusterProfiler")
    import rpy2.robjects as robjects
    go = topGO.runTest()
    
    robjects.r['']
    topGO.
    
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]

    #human importr(org.Hs.eg.db)


    clusterProfiler
    
    clusterProfiler, DOSE,
    

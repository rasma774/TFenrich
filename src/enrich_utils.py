import numpy as np
import scipy.stats as sts
import pandas as pd


__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2020 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

def _sortsets(db):
    # TODO: add tests that db is in ['KEGG', 'REACTOME', 'ALL']
    gene_lists = {}
    pw = __file__.split('/src')[0] 
    f = open(pw + '/data/gene_annotations/gene_annot/c2.all.v7.1.symbols.gmt')
    for line in f:
        line = line.replace('\n', '').split('\t')
        if (line[0][:len(db)] == db) or (db.upper() =='ALL'):
            genes = np.array(line[2:])
            gene_lists[line[0]] = genes
    f.close()
    return gene_lists

def _calc_fisher(gene_lists, genes_tmp, ngenes_thresh=10):
      
    # TODO: we really need to do some good selection for D!!!
    # Is this a good way?
    tmp = []
    [tmp.append(gene_lists[x]) for x in gene_lists.keys()]
    tmp.append(genes_tmp)
    flat_list = np.array([item for sublist in tmp for item in sublist])
    nunique = np.unique(flat_list).shape

    
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
        if len(gene_lists[hallmark]) < ngenes_thresh:
               continue
        A = np.in1d(gene_lists[hallmark], genes_tmp).sum()
        B = len(gene_lists[hallmark]) - A
        C = len(genes_tmp) - A
        D = nunique  - (A + B + C)

        OR[hallmark], p[hallmark] = sts.fisher_exact([[A, B], [C, D]], alternative='greater')
    res = pd.DataFrame([OR, p]).transpose()
    res.columns = ['OR', 'p']
    return res


def set_enrichments(gene_set, mult_test_corr=None, db='GO', FDR=0.05, ):
    """
    

    Parameters
    ----------
    gene_set : pd.array
        DESCRIPTION.
    db : str
        suggestions include, but are not limited to
        {'REACTOME', 'KEGG', 'BIOCARTA', 'PID', 'GO, 'ALL'}. The default is KEGG
    FDR : float, optional
        False discovey rate acc BenjaminiHochberg. 0 < FDR < 1. The default is 0.05.
    
    Returns
    -------
    enrichment analysis .

    """
    
    pw = __file__.split('/src')[0] 

    if type(db) is not str:
        gene_lists = db
    elif db.upper()  == 'GO':
        gene_lists = pd.read_pickle(pw + '/data/pickles/go_terms.p')
    elif db.upper() == 'GWAS':
        gene_lists = pd.read_pickle(pw + '/data/pickles/gwas.p')
    elif (db == 'KEGG') or (db == 'REACTOME'):
        # TODO: make this to a pickle too, and add 'ALL' as option for both
        gene_lists = _sortsets(db)
    else:
        raise ValueError('db not specified correctly, should be either dict, or string with values "GO", "GWAS", "KEGG", or "REACTOME"')
        
        
        
    res = _calc_fisher(gene_lists, gene_set)
    
    index_sort = np.argsort(res.p)
    res = res.iloc[index_sort, :]
    
    if not mult_test_corr is None:
            res['FDR'] = mult_test_corr(res.p, FDR=FDR)
    return res


# =============================================================================
#         
# def enrichr(genes):
#     import rpy2.robjects.packages as rpackages
#     from rpy2.robjects.packages import importr
#     #import rpy2
#     from rpy2 import robjects as robj
#     utils = importr('utils')
#     if not rpackages.isinstalled("org.Hs.eg.db"):
#         utils.install_packages("org.Hs.eg.db")
#         
#     if not rpackages.isinstalled("clusterProfiler"):
#         utils.install_packages("clusterProfiler")
#     importr("clusterProfiler")
#     import rpy2.robjects as robjects
#     go = topGO.runTest()
#     
#     robjects.r['']
#     topGO.
#     
#     names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
# 
#     #human importr(org.Hs.eg.db)
# 
# 
#     clusterProfiler
#     
#     clusterProfiler, DOSE,
#     
# =============================================================================


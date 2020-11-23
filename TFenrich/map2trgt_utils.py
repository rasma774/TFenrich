#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 13:34:37 2020

@author: rasmus
"""

import pandas as pd
import numpy as np

import stat_utils

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2020, Linköping'
__contact__ = 'rasma774@gmail.com'

def trrust_genes(TFs, weighted=False):
    # TODO: add something that corrects for genes that map to multiple TFs in input
    print('WARNING: should add something that has to do with number of times each target gene is used')
    """
    

    Parameters
    ----------
    TFs : list
        List of transcription factors to map to target genes.

    Returns
    -------
    List of TRRUST target genes.

    """
    # Load TRRUST
    print('loading TRRUST')
    TRRUST = pd.read_csv(
        '../data/TRRUST/trrust_rawdata.human.tsv', 
        sep='\t', 
        header=None,
        )
    print('Done')
    
    # As of now, we dont use the direction or publications
    TRRUST = TRRUST.iloc[:, :2]
    
    
    in_TFs = np.in1d(TRRUST[0].unique(), TFs)
    in_TRRUST = np.in1d(TFs, TRRUST[0].unique())


    print(str(100*np.sum(~in_TRRUST)/len(in_TRRUST)) + '% of TFs are not in TRRUST')
    print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of TRRUST TFs were in the TF list')
    
    target_genes = TRRUST[TRRUST[0].isin(TFs)][1].values
    unique_targets, counts = np.unique(target_genes, return_counts=True)
    
    if ~weighted:
        return unique_targets
    
    targets = pd.DataFrame((unique_targets, counts)).transpose()
    targets = targets.set_index(targets.columns[0]).iloc[:, 0]
    return targets


def correlation_genes(TFs, multiple_testing_correct=True, thresh=0.95):
    """
    

    Parameters
    ----------
    TFs : list
        List of transcription factors to map to target genes..

    Returns
    -------
    Panda series of correlating target genes summed over TFs.

    """
    print('loading corr')
    corr = pd.read_pickle('../data/pickles/correlations.p')
    corr = corr.set_index(corr.columns[0])
    print('Done')


    in_TFs = corr.index.isin(TFs)
    in_corr = np.in1d(TFs, corr.index)

    print(str(100*np.sum(~in_corr)/len(in_corr)) + '% of TFs are not found')
    print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of correlation table TFs were in the TF list')
  
    
    # Since the self-correlation is one, we need to remove the input TFs from 
    # the set
    corr_out = corr.iloc[:, ~corr.columns.isin(TFs)]
    target_genes = corr_out[in_TFs].abs().sum()

    if not multiple_testing_correct:
        return target_genes
    
    cval_dist = []
    for _ in range(100):
        randtfs = np.random.choice(corr.index, size=(in_corr.sum()), replace=False)
        ctmp = corr.iloc[:, ~corr.columns.isin(randtfs)][corr.index.isin(randtfs)].abs().sum()
        cval_dist.append(np.sort(ctmp)[int(len(ctmp)*thresh)])
    
    target_genes_adj = target_genes[target_genes >= np.median(cval_dist)]
    return target_genes_adj
    

def STRING_ppi(TFs, FDR=0.95, Npermut=100):
    """
    

    Parameters
    ----------
    TFs : TYPE
        DESCRIPTION.
    multiple_testing_correct : TYPE, optional
        DESCRIPTION. The default is True.
    thresh : TYPE, optional
        DESCRIPTION. The default is 0.95.

    Returns
    -------
    None.

    """
    
    print('loading STRING PPI...')
    ppi = pd.read_pickle('../data/pickles/string_links.p')
    print('Done')    
    
    in_TFs = np.in1d(ppi.index.unique(), TFs)
    in_STRINGdb = np.in1d(TFs, ppi.index.unique())

    print(str(100*np.sum(~in_STRINGdb)/len(in_STRINGdb)) + '% of TFs are not in STRINGdb')
    print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of STRINGdb TFs were in the TF list')
    
    summed_score = ppi[ppi.index.isin(TFs)].groupby('target_SYMBOL')['combined_score'].sum()
    summed_score = pd.DataFrame(summed_score)


    # TODO: The test here is not stringent at all. The more TFs, the more power,
    # and the more targets we get. Tested random 400 TFs, got 6500 significant genes.
    # But may be reasonable... ~25% of all possible TFs gave ~25% of all genes. 
    
    if FDR == -1:
        return summed_score
    else:    
        BH_sign = stat_utils._stringdb_bootstrap(summed_score.copy(), # to know what our gene scores are
                                                 ppi, # to get random gene scores
                                                 in_STRINGdb.sum(), # to know how many random TFs to draw
                                                 FDR=FDR, 
                                                 N=Npermut, # number of permutations
                                                 )
        return summed_score[BH_sign]
    
    
        
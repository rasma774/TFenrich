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

def trrust_genes(TFs, weighted=True):
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

    in_TFs = TRRUST.iloc[:, 0].isin(TFs)
    in_TRRUST = np.in1d(TFs, TRRUST[0].values)

    print(str(100*np.sum(~in_TRRUST)/len(in_TRRUST)) + '% of TFs are not in TRRUST')
    print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of TRRUST TFs were in the TF list')
    
    target_genes = TRRUST[in_TFs][1].values
    unique_targets, counts = np.unique(target_genes, return_counts=True)
    
    if ~weighted:
        return unique_targets
    
    targets = pd.DataFrame((unique_targets, counts)).transpose()
    targets = targets.set_index(targets.columns[0]).iloc[:, 0]
    return target_genes


def correlation_genes(TFs, multiple_testing_correct=True):
    # TODO :  add several cut-offs at the correlation level, no good way of
    # doing this
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
    corr = pd.read_csv(
        '../data/gene_correlations/correlations.csv', 
        index_col=0,
        )
    print('Done')


    in_TFs = corr.index.isin(TFs)
    in_corr = np.in1d(TFs, corr.index)

    print(str(100*np.sum(~in_corr)/len(in_corr)) + '% of TFs are not found')
    print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of TRRUST TFs were in the TF list')
  
    
    # Since the self-correlation is one, we need to remove the input TFs from the set
    corr = corr.iloc[:, ~corr.columns.isin(TFs)]
    
    target_genes = corr[in_TFs].abs().sum()

    return target_genes

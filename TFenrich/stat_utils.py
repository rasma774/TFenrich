#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 16:16:15 2020

@author: rasmus
"""

import numpy as np
import scipy.stats as sts

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2020, Linköping'
__contact__ = 'rasma774@gmail.com'

def benjaminihochberg_correction(p, FDR=0.05):
    """
    

    Parameters
    ----------
    p : np.array
        array of p-values of independent tests.
    FDR : float, optional
        False discovey rate. 0 < FDR < 1. The default is 0.05.

    Returns
    -------
    passes_FDR
        A vector of same length as input parameter p, with bool value True
        where the test passed a BH FDR correction.

    """
    sorted_p = np.sort(p)
    rank = np.arange(1, len(p)+1)

    BH_crit = (rank/rank[-1])*FDR

    if np.any(sorted_p < BH_crit):
        thresh = sorted_p[sorted_p < BH_crit][-1]
        return p <= thresh
    else:
        return p < 0
    

def _stringdb_bootstrap(summed_score, ppi, nTFs, FDR=0.05, N=100):
    unique_string_tfs = ppi.index.unique()
    for _ in range(N ):
        TFrand = np.random.choice(unique_string_tfs, size=nTFs, replace=False)
        random = ppi[ppi.index.isin(TFrand)].groupby('target_SYMBOL')['combined_score'].sum()
        summed_score[_] = random

    # Nan here just means not found, so should be 0
    summed_score[np.isnan(summed_score)] = 0
    
    means = summed_score.iloc[:, 1:].mean(1)
    stds = summed_score.iloc[:, 1:].std(1)
    Z = (means - summed_score.combined_score)/stds
    p = sts.norm.cdf(Z)
    is_sign = benjaminihochberg_correction(p, FDR=FDR)
    return is_sign

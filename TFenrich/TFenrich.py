#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:30:20 2020

@author: rasmus

# TODO: add docstrings when we see that everything works

"""

import enrich_utils
import map2trgt_utils
import stat_utils

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2020, Linköping'
__contact__ = 'rasma774@gmail.com'

# TODO: we need some way to handle rankings, and not just boolean in/not in list
#       in the enrichment calculations.
# TODO: we need to find some good way to get less genes from STRINGdb 
# TODO: if the program is run at the first time, assemble the correlations table
# TODO: add silent option
class TFenrich:
    def __init__(self, 
                 TFs, 
                 mapmethod='corr', 
                 multiple_testing_correction='BenjaminiHochberg',
                 silent=False):
        """
        # TODO: add description here        

        Parameters
        ----------
        TFs : list or dict
            If list, a list of the TFs that are subject to the analysis.
            If a dict, a dict of TFs and weights. Only applicable if mapmethod
            is 'deep'
        mapmethod : {'deep', 'corr', 'TRRUST'}, optional
            Select the type of mapping from TFs to targets. The types are 
            either of the classes
            deep: Uses a deep learning model that maps TF expression to 
                targets 
            corr: Uses a matrix of gene expression correlation
            TRRUST: Uses the TRRUST database of known physical interactions
            all: All of the available mappings
            The default is 'deep'.

        Returns
        -------
        None.

        """
        
        self.TFs = TFs
        self.silent = silent
        if multiple_testing_correction == 'BenjaminiHochberg': # use this unless otherwise told
            self.multtest_fun = stat_utils.benjaminihochberg_correction
        
        self.mapmethod = mapmethod
        if self.mapmethod == 'corr':
            self.target_genes = map2trgt_utils.correlation_genes(TFs, silent).index
        elif self.mapmethod == 'TRRUST':
            self.target_genes = map2trgt_utils.trrust_genes(TFs, silent)
        elif self.mapmethod == 'PPI':
            self.target_genes = map2trgt_utils.STRING_ppi(TFs, silent)
        else:
            raise ValueError('mapmethod \'' + self.mapmethod + '\' not defined')
            
            
    def downstream_enrich(self, mult_test_corr='same', db='GO', FDR=0.05):

        # If not specified use global
        if mult_test_corr == 'same':
            mult_test_corr = self.multtest_fun
            
        res = enrich_utils.set_enrichments(self.target_genes,
                                           mult_test_corr=mult_test_corr,
                                           db=db, 
                                           FDR=FDR,
                                           )
        self.pathway_enrichments = res
    
    # TODO: Add some plotting function here!
    #   Suggestion: dotplot?
    
        
    
        
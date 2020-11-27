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
import plot_utils

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Rasmus Magnusson, 2020, Linköping'
__contact__ = 'rasma774@gmail.com'

# TODO: we need some way to handle rankings, and not just boolean in/not in list
#       in the enrichment calculations.
# TODO: we need to find some good way to get less genes from STRINGdb 
# TODO: if the program is run at the first time, assemble the correlations table
# TOOD: add argparse
# TODO: Handle p-values as negative log10, such that the Fisher test can be estimated
# TODO: add the DisGenet database to compare with

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
            self.target_genes = map2trgt_utils.correlation_genes(TFs, silent=silent).index
        elif self.mapmethod == 'TRRUST':
            self.target_genes = map2trgt_utils.trrust_genes(TFs, silent=silent)
        elif self.mapmethod == 'PPI':
            self.target_genes = map2trgt_utils.STRING_ppi(TFs, silent=silent)
        else:
            raise ValueError('mapmethod \'' + self.mapmethod + '\' not defined')
            
            
    def downstream_enrich(self, mult_test_corr='same', db='GO', FDR=0.05):
        """
        

        Parameters
        ----------
        mult_test_corr : TYPE, optional
            DESCRIPTION. The default is 'same'.
        db : TYPE, optional
            DESCRIPTION. The default is 'GO'.
        FDR : TYPE, optional
            DESCRIPTION. The default is 0.05.

        Returns
        -------
        None.

        """

        # If not specified use global
        if mult_test_corr == 'same':
            mult_test_corr = self.multtest_fun
            
        res = enrich_utils.set_enrichments(self.target_genes,
                                           mult_test_corr=mult_test_corr,
                                           db=db, 
                                           FDR=FDR,
                                           )
        self.enrichments = res
    
    def plot(self, 
             savename=None,
             plot_Ntop='all',
             textlength=None,
             sorton='OR',
             remove_non_FDR=True,
             ):
        """
        

        Parameters
        ----------
        savename : TYPE, optional
            DESCRIPTION. The default is None.
        number_to_plot : TYPE, optional
            DESCRIPTION. The default is 'all'.

        Returns 
        -------
        None.

        """
        
        if plot_Ntop == 'all':
            plot_Ntop = self.enrichments.shape[0]

        # TODO: add kwargs?
        f, ax = plot_utils.plot_pvals(self.enrichments.copy(), 
                                      savename, 
                                      plot_Ntop,
                                      textlength=textlength,
                                      sorton=sorton,
                                      remove_non_FDR=remove_non_FDR,
                                      )
        return f, ax
    
        



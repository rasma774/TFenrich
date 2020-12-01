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
#       in the enrichment calculations. Would, howvever, not work with 
#       the binding sites such as TRRUST. OTherwise use KS test
# TODO: if the program is run at the first time, assemble the correlations table
# TOOD: add argparse
# TODO: Handle p-values as negative log10, such that the Fisher test can be estimated
# TODO: add the DisGenet database to compare with
# TODO: Do an anlysis on the clutering of random TFs and the three downstream mappings .
# TODO: get a list of all citations in a dict, and append it to the method. Put 
#       all citations in a textfile
# TODO: I think i should remove TRRUST, there are too many biases atm...
# TODO: take decision on whether to only have the correlation-based analysis.
#       Seems like the PPi is heavily biased.

class TFenrich:
    def __init__(self, 
                 TFs, 
                 mapmethod='corr', 
                 silent=False,
                 top_n_genes=None):
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
        assert len(TFs) > 0
        
        self.TFs = TFs
        self.silent = silent

        if mapmethod == 'corr':
            self.mapmethod = map2trgt_utils.correlation_genes
            
        self.target_genes = self.mapmethod(TFs, 
                                           silent=silent,
                                           top_n_genes=top_n_genes
                                           ).index.values
            
            
    def downstream_enrich(self, 
                          db='GO', 
                          FDR=0.05,
                          multiple_testing_correction='BenjaminiHochberg',
                          ):
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

                
        if multiple_testing_correction == 'BenjaminiHochberg': # use this unless otherwise told
            self.multtest_fun = stat_utils.benjaminihochberg_correction
        else:
            # Here, we let the user define the testing function
            self.multtest_fun = multiple_testing_correction

        res = enrich_utils.set_enrichments(self.target_genes,
                                           mult_test_corr=self.multtest_fun,
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
        f, ax = plot_utils.plot_res(self.enrichments.copy(), 
                                      savename, 
                                      plot_Ntop,
                                      textlength=textlength,
                                      sorton=sorton,
                                      remove_non_FDR=remove_non_FDR,
                                      )
        return f, ax
    
        



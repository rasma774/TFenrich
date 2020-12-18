from src import enrich_utils
from src import map2trgt_utils
from src import stat_utils
from src import plot_utils
from src import parse_utils
from src import citation_handler

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2020 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'
__LICENSE__ = 'GNU Affero General Public License v3.0'
__version__ = '0.01'
__cite__ = ''


# TODO: if the program is run at the first time, assemble the correlations table
# TODO: Handle p-values as negative log10, such that the Fisher test can be estimated
# TODO: add the DisGenet database to compare with

class TFenricher:
    def __init__(self, 
                 TFs, 
                 mapmethod='corr', 
                 silent=False,
                 top_n_genes=None):
        """

        Parameters
        ----------
        TFs : list or dict
            If list, a list of the TFs that are subject to the analysis.
            If a dict, a dict of TFs and weights. Only applicable if mapmethod
            is 'deep'
        mapmethod : function, optional
           

        Returns
        -------
        None.

        """
        assert len(TFs) > 0
        
        self.TFs = TFs.copy()
        self.silent = silent
        
        # For the references
        self.used_methods = []
                    

        if mapmethod == 'corr':
            self.mapmethod = map2trgt_utils.correlation_genes
            self.used_methods.append('corrs')
        else:
            self.mapmethod = mapmethod
    
        self.target_genes = self.mapmethod(TFs, 
                                           silent=silent,
                                           top_n_genes=top_n_genes
                                           )
            
            
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

        if type(db) == str:
            self.used_methods.append(db)
        
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
             padding=0.7,
             tick_font_size=14):
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

        # TODO: add kwargs
        f, ax = plot_utils.plot_res(self.enrichments.copy(), 
                                      savename, 
                                      plot_Ntop,
                                      textlength=textlength,
                                      sorton=sorton,
                                      remove_non_FDR=remove_non_FDR,
                                      padding=padding,
                                      tick_font_size=tick_font_size,
                                      )
        return f, ax
    
    def cite(self):
        """
        Write the references to the databases used, and the TFenricher reference

        Returns
        -------
        None.

        """
        citation_handler.write_citations(self.used_methods)
        
        
        
if __name__ == '__main__':
    args = parse_utils.parse()
    # Map TFs to targets
    enr = TFenricher(args.tfs, silent=args.silent, top_n_genes=args.ngenes)
    
    # Calculate the overlaps between putative downstream genes and gene sets 
    enr.downstream_enrich(db=args.db, FDR=args.FDR)
    
    if args.plotname != '-1':
        enr.plot(savename=args.plotname, plot_Ntop=args.plot_n_top[0])
    
    enr.enrichments.to_csv(args.results_savename)
        

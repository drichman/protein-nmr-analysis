Example code and data set for analysis of hydrogen-deuterium exchange 
in proteins. By Dan Richman, 2015 March 31.


USAGE

Run from the Unix command line with "python hdx_analysis.py" as long 
as your python is Python 3. (The other file, nmrfn.py, is a module 
containing some functions called by hdx_analysis.py.)

There will be a few seconds of seeming inactivity while modules load,
then excluded residues will be printed to the screen, then output files 
will be written.


EXPECTED OUTPUT FILES

hdx_rates.txt - columns:
    (1st col) : pH
    (2nd col) : residue number
    dec_frac  : signal decay ratio (end height over start height)
    k_obs     : observed exchange rate (1/s)
    k_error   : observed exchange rate error
    er_frac   : observed exchange rate fractional error
    dg        : delta G derived from the exchange rate (kcal/mol)
    dger      : delta G error

hdx_dg.txt : derived delta G for each residue at each of 3 pH conditions

hdx_dg_bins.txt : residues sorted by which range their delta G values
                  fall into (eg 0-3 kcal/mol, 6-9 kcal/mol) at each pH

hdx_slopevals.txt : delta G change as pH drops from pH* 4.5 to 3.9
                    (deltaG @ pH* 4.5 minus deltaG @ pH* 3.9 delta)/(4.5-3.9)

hdx_dg_slope_bins.txt : residues sorted by how big the delta G change is

hdx_dg_slopes_plot.pdf : delta G slopes plotted vs sequence, annotated with 
                         locations of surface Asp and Glu residues in SNase 
                         and with the approximate change in global delta G 
                         between pH 4.5 and 3.9



ABOUT THE INPUT DATA

Peak heights tables from sparky, one for each of 3 pH conditions 
(actually pH* ie uncorrected meter reading). Row labels are amino acid residue, 
column labels are datetime stamps.
(The protein is Delta+PHS staph nuclease.)

Intrinsic exchange rates simulated by the SPHERE program online.


ABOUT THE CODE

Extracts a hydrogen exchange rate for each residue by fitting 
decaying exponentials to the peak height vs datetime data. Automates reads 
the pH condition from the peak height file name and contains functions for 
deriving local stability from the exchange rate, comparing the results 
between pH conditions, and writing tables and plots.


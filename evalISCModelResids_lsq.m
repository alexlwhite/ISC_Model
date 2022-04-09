function resid = evalISCModelResids_lsq(params)

global dataToFit; 

fitdata = ISCModel(params);

resid = dataToFit-fitdata;

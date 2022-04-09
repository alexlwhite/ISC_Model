function sqre = evalISCModelResids_fmin(params)

global dataToFit; 

fitdata = ISCModel(params);
resid = dataToFit(:)-fitdata(:);
sqre = sum(resid(:).^2);
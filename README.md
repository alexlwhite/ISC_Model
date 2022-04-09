# ISC_Model
Data and Matlab code for the "Independent Systems with Competition" (ISC) model from White, Rolfs and Carrasco (2015)
Running the script RunISCModel will load in real data, and generate predicted d' levels using the function ISCModel. It saves the predictions in a mat file ISCModelFit_{date} and makes a figure DataWithISCModelPreds.eps. 

The real psychophysical data are contained in the Matlab file FBAxSADPrimes.mat
matrices X1Ds and X2Ds, of d' accuracy measures.
% Each of those is a 3x3xN matrix, where N is the number of subjects. 
% As described in the struct "dimensionLabels" (also in the .mat file),
% dimension 1 is the spatial cue validity, dimension 2 is feature cue
% validity, and dimension 3 is for individual subjects. 

The function ISCModel takes in 6 parameters and ouputs predicted d' data. 

The functions evalISCModelResids_fmin and evalISCModelResids_lsq are for actually running optization routines (either fminsearch) 
or lsqnonlin to find best-fitting parameters. 

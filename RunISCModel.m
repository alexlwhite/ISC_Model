%%%%% "Independent Systems with Competition" (ISC) Model for the Interaction of Spatial and Feature-based Attention %% %%
% % % by Alex White
% % % January, 2014 (cleaned up April 2022)
% 
% This fits the ISC model to psychophysical data, as reported in: 
%﻿White, A. L., Rolfs, M., & Carrasco, M. (2015). Stimulus competition mediates the joint effects of spatial and feature-based attention. Journal of Vision, 15, 1–21.
%     https://doi.org/10.1167/15.14.7
% 
% Note: this script may no longer work to efficiently run an optimization
% routine on the parameters (but it tries if runFit==1), and it does not include code to run a grid
% search through parameters as we actually did to find the best-fitting
% parameters. As currently set up, this script assumes those parameters and
% just computes and plots the model predictions. 
%% Set up parameters
clear; close all; home;

global ntrials
global dataToFit
global globalMean
global attnIs

%ntrials: how many trials are simulated
ntrials = 500000;

%runFit: 1 = truly fit; 0 =  just test, by computing the predictions of one set of params;
runFit = 0;

optFun = 2; %Optimization function: 1=lsqnonlin; 2=fminsearch

%The attnIs matrix has one row for each attention condition, and each
%item in that row is the index number for the corresponding position in
%the final output

attnIs(1,:) = [5 5 5 5]; %neutral spatial, neutral feature
attnIs(2,:) = [2 2 8 8]; %neutral spatial, invalid and valid feature
attnIs(3,:) = [4 6 4 6]; %invalid and valid spatial, with neutral feautre
attnIs(4,:) = [1 3 7 9]; %invalid spatial, invalid feature; valid spatial, invalid feature; invalid spatial, valid feature; valid spatial, valid feature; 

%Here's another way to describe the indices in predicted d-prime levels, 
% which correspond to values of attnIs. 
% The dprime matrices, dataToFit model predictedDs, have one row for each experiment, with each item
%ordered as follows:
%   1     2     3     4     5     6     7     8     9
%[FI_SI FI_SN FI_SV FN_SI FN_SN FN_SV FV_SI FV_SN FV_SV]
% where F=feature-based, S=spatial, and I=invalid,N=neutral,V=valid.

%% Load and plot data

%load a mat file that has individual subject actual data from both
%experiments: matrices X1Ds and X2Ds, of d' accuracy measures.
% Each of those is a 3x3xN matrix, where N is the number of subjects. 
% As described in the struct "dimensionLabels" (also in the .mat file),
% dimension 1 is the spatial cue validity, dimension 2 is feature cue
% validity, and dimension 3 is for individual subjects. 
load('FBAxSADPrimes.mat');

nCueConds = size(X1Ds,1);
nDs = nCueConds^2; %number of data points modeled

close all;

fig = figure;

dataToFit = NaN(2, nDs);
for xi=1:2
    
    if xi==1
        %pull out individual 
        ds = X1Ds;
        ttl = 'Expt 1: With distractors';
        addLabels = true;
    else
        ds = X2Ds;
        ttl = 'Expt 2: No distractors';
        addLabels = false;
        
    end
    
    mD = mean(ds,3); %across-subject means
    
    %compute normalized within-subject SEMs
    smeans = mean(mean(ds,1),2);
    smeansrep = repmat(smeans,[nCueConds nCueConds 1]);
    dsnorm = ds-smeansrep;
    SEMs = nDs/(nDs-1)*std(dsnorm,0,3)/sqrt(size(ds,3));
    
    subplot(1,2,xi); hold on;
    
    barCtrs = plotFBAxSAData(mD, SEMs, addLabels);
    title(ttl);
    
    dataToFit(xi, :) = mD(:)';
end

%global mean of those data, for rescaling model preds
globalMean = mean(dataToFit(:));

%% Fit

paramLabs = {'base','k','sig_s','sig_f','SNR1','SNR2'}; %free parameters
nParams = numel(paramLabs);

startParams = zeros(1,nParams);
startParams(1) = 1.0148;   % response to baseline saturation
startParams(2) = 6.394;    % k, the exponent
startParams(3) = 0.0294;   % sig_s, spatial cueing addition for dots with saturation change
startParams(4) = 0.1218;   % sig_f, feature cueing addition for dots with saturation change
startParams(5) = 6.931;    % r1, the starting SNR in neutral-neutral condition of X1
startParams(6) = 1.025;    % r2, the starting SNR in neutral-neutral condition of X2

%bounds on parameters:
maxattn = 0.9;
lb = [1   3   0          0    0.1  0.1];
ub = [25  15  maxattn  maxattn 20   20];


if runFit %actually fit 
    if optFun == 1
        options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-12,'TolX',1e-13,'DiffMinChange',1e-12);
        [fitParams,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]= lsqnonlin(@evalISCModelResids_lsq,startParams,lb,ub,options);
    else
        options = optimset('Display','iter','MaxFunEvals',400,'MaxIter',500,'DiffMinChange',1e-6);%,'TolFun',1e-12,'TolX',1e-13);
        startParams1 = startParams;
        for fitrun = 1:5
            [fitParams,FVAL,EXITFLAG,OUTPUT] = fminsearch(@evalISCModelResids_fmin,startParams,options);
            fprintf(1,'\nfminsearch run %i, startParams =', fitrun);
            disp(startParams);
            fprintf(1,'\t\tfitParams  =  ');
            disp(fitParams);
            startParams = fitParams;
        end
    end
else %just test the starting parameters 
    fitParams = startParams;
end

%% Get the model predictions
[predictedDs, scaleFactor] = ISCModel(fitParams);

%computer r-squared 
resid = dataToFit(:)-predictedDs(:);
SSTot = sum((dataToFit(:)-mean(dataToFit(:))).^2);
SSE   = sum(resid.^2);
rSqr  = 1-SSE/SSTot;

%or, in simpler terms:
%rSqr2 = 1 - ( norm(dataToFit(:)-predictedDs(:)) / norm(dataToFit(:)-mean(dataToFit(:))) )^2

%% plot the model predictions on top of the daa
for xi=1:2
    dmx = predictedDs(xi,:);
    
    figure(fig); subplot(1,2,xi); hold on;
    set(get(gca,'Legend'),'AutoUpdate','off');

    %somewhat silly way of plotting each point as a bullseye
    plot(barCtrs,dmx(:),'k.','MarkerSize',28);
    plot(barCtrs,dmx(:),'w.','MarkerSize',24);
    plot(barCtrs,dmx(:),'k.','MarkerSize',12);
    
end

%add some text
textSz = 11;
dylims = [0 2.05]; %y-axis limits 

height=1;
hdiff=0.06;

%print params
for pi = 1:nParams
    height=height-hdiff;
    text(0.8, dylims(2)*height,sprintf('%s=%.2f',paramLabs{pi},fitParams(pi)),'FontSize',textSz);
end
height = height-hdiff;
text(0.8, dylims(2)*height,sprintf('%s=%.2f','scaleF',scaleFactor),'FontSize',textSz);

height=height-hdiff*2;
text(9.6, dylims(2)*0.94,sprintf('r^2=%.2f',rSqr),'FontSize',textSz);


for spi=1:2
    subplot(1,2,spi); box off; axis square;
end

set(gcf,'Position',[5 5 700 500]);
set(gca,'FontSize',18);

if runFit
    disp(OUTPUT);
end


fName = 'DataWithISCModelPreds.eps';
exportfig(gcf,fName,'Format','eps','bounds','tight','color','rgb','LockAxes',0,'FontMode','scaled','FontSize',1);

save(sprintf('ISCModelFit_%s.mat',date),'fitParams','predictedDs','rSqr')


%% function [dmNorm, scaleFactor, dm] = ISCModel(params)
% This function predicts d-prime accuracy levels given some parameters for
% the Independent Systems with Competition (ISC) model, as described in 
% White, Rolfs & Carrasco, JOV 2015. (https://doi.org/10.1167/15.14.7) 
% 
% Inputs: 
% - params: a 1x6 vector of parameters. 
%    1 = response magnitude to dots of baseline saturation level, u_base
%    2 = exponent, k
%    3 = spatial attention modulation, S
%    4 = feature attention modulation, F 
%    5 = incremental response to saturation increment in Expt 1 (delta_1) 
%    6 = incremental response to saturation increment in Expt 2 (delta_2) 
% 
% Outputs: 
% - dsScaled: predicted dprime levels, after multiplying raw predictions (dsNormed) by same scaleFactor
%   such that the mean of dsScales is the same as mean of true data points. 
% - scaleFactor: the number by which raw predictions (dsNormed) got multipled 
% - dsNormed: raw predictions before rescaling. 

function [dsScaled, scaleFactor, dsNormed] = ISCModel(params)


global globalMean
global ntrials
global attnIs

b = params(1); %mean baseline response to dots with no increment 
k = params(2); %exponent for "competition"
s = params(3); %spatial attention modulation 
f = params(4); %feature attention modulation 
r1 = params(5);%sensory increment response in X1 
r2 = params(6);%sensory increment response in X2

incrementNeutralResponse1 = b + r1; %response to incremented dots in Expt 1 (four increments at a time)
incrementNeutralResponse2 = b + r2; %response to incremented dots in Expt 2  only one increment at a time)


%Four attentional conditions of feature * space attention: 
% neutral*neutral, cued*neutral, neutral*cued, cued*cued

attnGains(1,:) = [1 1 1 1]; %neutral spatial, neutral feature.
attnGains(2,:) = 1+[-f -f f f ]; %spatial neutral, feature invalid or valid 
attnGains(3,:) = 1+[-s s -s s]; %feature neutral, spatial invalid or valid
attnGains(4,:) = attnGains(2,:)+attnGains(3,:)-1;  %the cross-combos of feature and spatial each invalid or valid. Need to subtract 1 so each item is just (1+f+s)

%these can't go below 0, because you cant have negative means
attnGains(attnGains<0)=0;

%Indices: the output has one row for each experiment, with each item
%ordered as follows(F=feature-based, S=spatial, I=invalid,N=neutral,V=valid):
%   1     2     3     4     5     6     7     8     9
%[FI_SI FI_SN FI_SV FN_SI FN_SN FN_SV FV_SI FV_SN FV_SV]

%So each of the four attentional conditions listed above will be used to
%measure a few of those conditions.

nCond = size(attnIs, 1);
nDs  = 9; %number of individual data points

%predicted  dprime levels, "normed" before scaling
dsNormed = zeros(2,nDs);

%baseline mean response to all dot fields:
mus0 = ones(ntrials,4,2)*b; %all start with baseline level b

for xi = 1:2
    %For Expt 2 (four saturation changes at once), set the means of saturation responses.
    %(Expt 3 is done later, for 1 dot field at a time). 
    if xi==1 
        % Assume that the increments are all in the upper half. 
        %That didn't happen in real experiment, but it
        %doesn't matter for the model.
        mus1 = mus0;
        mus1(:,:,2) = incrementNeutralResponse1;
    end
    
    for attncond = 1:nCond %for each attention condition
        
        ais = attnIs(attncond,:); %which points we have in this attention condition
        
        for ai = unique(ais)
            
            %% Stage 1: Encoding + Attentional modulation.
            % First, set the means of the distributions of saturation estimates
            % for each sub-field of dots (some of which may have a saturation increment)
            
            %si is the index of item with the saturation change / the item that is post-cued
            si = find(ais==ai);
            si = si(1);
            
            %set stimulus response means. For X1, that's done already as
            %mus1. For X2, we need to set just 1 to have the increment: 
            if xi==2
                %just 1 gets the increment
                mus1 = mus0;
                mus1(:,si,2) = incrementNeutralResponse2; %one gets the increment
            end
            
            %Add attention effects: 
            % multiply all the mean saturation responses by attentional gain changes:
            mus = mus1.*repmat(attnGains(attncond,:),[ntrials 1 2]);
            
            %Now actually draw saturation responses, called "xs' from Gaussian
            %distributions with those means, and SD 1.
            xs = randn(ntrials,4,2) + mus; %assume Std Dev of noise is 1
            %size of xs: ntrials x 4 x 2
            
            %don't allow negative responses
            xs(xs<0) = 0;
            
            %% Stage 2: Normalization
            expds = xs.^k;
            xsn = expds./repmat(sum(sum(expds,2),3),[1 4 2]);
            
            %% Stage 3: Difference between top and bottom dots saturation estimate:
            ds = diff(xsn,1,3);
            
            %% Stage 4: Late noise
            ys = ds+randn(ntrials,4);
            
            %% Stage 4: Subject's explicit decision, which can be incorrect or correct: 
            corrects = squeeze(ys(:,si)>0);
            
            %avoid 100% accuracy by setting 1 response to be incorrect
            if all(corrects)==1, corrects(end)=0; end
            
            %% Compute dprime, assuming neutral criterion
            dsNormed(xi,ai) = 2*norminv(mean(corrects));
            
            
        end
    end
end

dsNormed(dsNormed<0)=0; %dprimes can't be below 0

%% Lastly, scale all points together to match mean of real data 
scaleFactor = globalMean/mean(dsNormed(:));
dsScaled = dsNormed*scaleFactor;



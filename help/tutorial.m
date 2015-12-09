%% Load the raw data
rawData = load('exampleData.mat'); % run tutorial_exampleData to generate this
nTrials = rawData.nTrials; % number of trials
trial=rawData.trial;
unitOfTime = 'ms';
binSize = 1; % TODO some continuous observations might need up/down-sampling if binSize is not 1!?

%% Specify the fields to load
n=neuroGLM('ms', binSize, rawData.param);

%% Convert the raw data into the experiment structure

n.addCovariateSpiketrain(trial, 'hist', 'sptrain', 'History Filter')
n.addCovariateSpiketrain(trial, 'coupling', 'sptrain2', 'Coupling from neuron 2')

%% Build 'designSpec' which specifies how to generate the design matrix
% Each covariate to include in the model and analysis is specified.


bs = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 10, n.binfun);
bs.normalize

%% Instantaneous Raw Signal without basis
n.addCovariateRaw(trial, 'LFP', 'LFP', bs);

%% Acausal Timing Event
bs = basisFactory.makeSmoothTemporalBasis('boxcar', 300, 8, n.binfun);
offset = -200;
n.addCovariateTiming(trial, 'saccade', 'saccade', [], bs, offset);

% %% Coherence
% % a box car that depends on the coh value
% bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 200, 10, n.binfun);
% bs.normalize
% stimHandle = @(trial) trial.coh * basisFactory.boxcarStim(binfun(trial.dotson), binfun(trial.dotsoff), binfun(trial.duration));
% 
% dspec = buildGLM.addCovariate(dspec, 'cohKer', 'coh-dep dots stimulus', stimHandle, bs);

% %% 2-D eye position
% bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 40, 4, binfun);
% dspec = buildGLM.addCovariateRaw(dspec, 'eyepos', [], bs);

%buildGLM.summarizeDesignSpec(dspec); % print out the current configuration

%% Compile the data into 'DesignMatrix' structure
trialIndices = 1:(nTrials-1); % use all trials except the last one
n.compileDesignMatrix(trial, trialIndices);


%% Get the spike trains back to regress against
y = n.getBinnedSpikeTrain(trial, 'sptrain', n.dm.trialIndices);

%% Do some processing on the design matrix
n.removeConstantCols;
colIndices = n.getDesignMatrixColIndices('LFP');
n.zscoreDesignMatrix([colIndices{:}]);

n.addBiasColumn; % DO NOT ADD THE BIAS TERM IF USING GLMFIT

%% Least squares for initialization
tic
wInit = n.dm.X \ y;
toc

%% Use matRegress for Poisson regression
% it requires `fminunc` from MATLAB's optimization toolbox


fnlin = @expfun; % inverse link function (a.k.a. nonlinearity)
lfunc = @(w)neglogli_tspline_poiss(w, n.dm.X, y, fnlin, 1); % cost/loss function

opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');

[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
wvar = diag(inv(hessian));

%% Alternative maximum likelihood Poisson estimation using glmfit
% [w, dev, stats] = glmfit(dm.X, y, 'poisson', 'link', 'log');
% wvar = stats.se.^2;

%% Visualize
ws = n.combineWeights(wml);
wvar = n.combineWeights(wvar);

fig = figure(2913); clf;
nCovar = numel(n.covar);
for kCov = 1:nCovar
    label = n.covar(kCov).label;
    subplot(nCovar, 1, kCov);
    plot(ws.(label).tr, ws.(label).data, ...
        ws.(label).tr, ws.(label).data+sqrt(wvar.(label).data), '--', ...
        ws.(label).tr, ws.(label).data-sqrt(wvar.(label).data), '--')
    title(label);
end

return

%{
%% Specify the model
hasBias = true;
model = buildGLM.buildModel(dspec, 'Poisson', 'exp', hasBias);

%% Do regression
[w, stats] = fitGLM(model, dm, y);
%}

%% Visualize fit
visualizeFit(w, model, dspec, vparam(1)); % ???

%% Simulate from model for test data
testTrialIndices = nTrial; % test it on the last trial
dmTest = compileSparseDesignMatrix(expt, dspec, testTrialIndices);

yPred = generatePrediction(w, model, dmTest);
ySamp = simulateModel(w, model, dmTest);

%% Validate model
gof = goodnessOfFit(w, stats, model, dmTest);
visualizeGoodnessOfFit(gof);

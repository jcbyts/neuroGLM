%% Generate example (raw) dataset for the tutorial

% REQUIRES stat toolbox for poissnrnd

nTrials = 500; % total number of trials

% preallocate structure in memory
trial = struct();
trial(nTrials).duration = 0; % preallocate

lambda = 15; % (baseline) firing rate per bin

% build glm kernels
maxLag = 1e3;
stimulusLags = (0:maxLag)'; % long (for LIP)

% make kernels an alpha function
fun = @(t, tau1, tau2) exp(-t/tau1).*(1 - exp(-t/tau2));
% fun = @(t, lag1, int1, mag2, lag2, int2) exp(- (t - lag1).^2 / (2*int1^2) ) - mag2 * exp( - (t - lag2).^2 / (2* int2^2));

kCon = -.01*fun(stimulusLags, 100, 20);
kDir = fun(stimulusLags, 100, 20);
sacOffset = 400;
kSac = stimulusLags/maxLag*15;
kSac(sacOffset:end) = 0;
figure(1); clf
plot(stimulusLags, kSac)
%%

for kTrial = 1:nTrials
    %%
    duration = 1550 + ceil(rand * 500);
    trial(kTrial).duration  = duration;
    % dots come on at a random time
    trial(kTrial).motionon    = 300 + ceil(rand * 500);
    % exponential random variable for motion duration
    trial(kTrial).motionoff   = trial(kTrial).motionon + 100 + ceil(min(exprnd(50), 200));
    motionDuration = trial(kTrial).motionoff-trial(kTrial).motionon+1;
    trial(kTrial).contrast    = sparse(trial(kTrial).motionon:trial(kTrial).motionoff, ones(motionDuration, 1), ones(motionDuration, 1), duration, 1);
    
    % instantaneous direction
    trial(kTrial).coh         = sign(rand - 0.5) * 2^ceil(rand*6)/100;
    % boxcar of height coh
    direction = trial(kTrial).coh * sparse(trial(kTrial).motionon:trial(kTrial).motionoff, ones(motionDuration, 1), ones(motionDuration, 1), duration, 1);
    % coherence "noise"
    motionNoise = randn(motionDuration, 1) * 5;
    direction(trial(kTrial).motionon:trial(kTrial).motionoff) = direction(trial(kTrial).motionon:trial(kTrial).motionoff) + motionNoise;
    trial(kTrial).direction = direction;
    
    % some sort of psychophysics
    d = sum(trial(kTrial).direction) + randn*2; 
    
    trial(kTrial).saccade   = trial(kTrial).motionoff - 50 + ceil(rand * sacOffset);
    
    trial(kTrial).choice    = sign(d);
    
    % build a neuron that responds to these stimuli in a stereotyped way
    xSac = sparse(ceil(trial(kTrial).saccade-sacOffset), 1, 1, duration, 1)*trial(kTrial).choice;
    lamRate = lambda + filter(kDir, 1, full(direction)) + filter(kCon, 1, full(trial(kTrial).contrast)) + filter(kSac, 1, full(xSac));
%     lamRate = exp(lamRate);
    figure(1); clf
    plot(lamRate)
    hold on
    plot(trial(kTrial).saccade*[1 1], ylim)
    plot(xSac)
    plot(direction)
    plot(trial(kTrial).contrast)
    title(trial(kTrial).coh)
    drawnow

    spks = rand(1, duration)<(lamRate/1e3);
    plot(spks*100, 'k')
    
    sptimes = find(spks);
    trial(kTrial).sptrain = sptimes;
    
    
%     trial(kTrial).sptrain2 = sort(rand(poissrnd(0.1 * lambda * duration), 1) * duration);
% 
%     trial(kTrial).sptrain = sort([trial(kTrial).sptrain; trial(kTrial).sptrain2 + 2]);
%     trial(kTrial).sptrain(trial(kTrial).sptrain > trial(kTrial).duration) = [];
end

param.samplingFreq = 1; % kHz
param.monkey = 'F99';

save('exampleData.mat', 'nTrials', 'trial', 'param')

%%



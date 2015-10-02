function model = doRegressionPoisson(X, Y, dspec, ndx, dt)
% Helper for poisson regression fits
% model = doRegressionPoisson(X, Y, dspec, k0, ndx)
% INPUTS:
%   X [nsamples x ndims] Stimulus
%   Y [nsamples x 1]     Count
%   dspec [struct] or [neuroGLM object]
%       .covar - array of covariates
%           .label
%           .edim
%           .desc
%       .edim  - total dimensionality
%       .unitoftime
%       .binSize
%       .model
%           .regressionMode
%           .bilinearMode
%           .bilinearRank
%           .bilinearCovariate


model=[];

k0=(X'*X + eye(size(X,2)))\(X'*Y);

assert(isfield(dspec, 'model'), 'model is a required field of dspec')
if ~isfield(dspec.model, 'regressionMode')
    dspec.model.regressionMode='ML';
end

if ~exist('ndx', 'var') || isempty(ndx)
    ndx = 1:numel(Y);
end

if ~exist('dt', 'var')
    dt = 1;
end

if isfield(dspec.model, 'optimOpts') && isa(dspec.model.optimOpts, optim.options.Fminunc)
    optimOpts=dspec.model.optimOpts;
else
    optimOpts = optimoptions(@fminunc, 'Display', 'iter', 'Algorithm','trust-region',...
        'GradObj','on','Hessian','on');
end

if isfield(dspec.model, 'nlfun') && isa(dspec.model.nlfun, 'function_handle')
    nlfun=dspec.model.nlfun;
else
    nlfun = @expfun; % canonical link @logexp1;
end

switch dspec.model.regressionMode
    case {'MLEXP', 'MLSR', 'MLE', 'ML'}
        lfun = @(w) neglogli_poiss(w, X(ndx,:), Y(ndx), nlfun, dt);
        [wmle, ~, ~, ~, ~, H] = fminunc(lfun, k0, optimOpts);
        model.khat      = wmle;
        model.fnlin     = nlfun;
        model.SDebars   = sqrt(diag(inv(H)));
    case {'Ridge', 'RIDGE'}
        [wRidge,rho,SDebars] = autoRegress_PoissonRidge(X(ndx,:),Y(ndx),nlfun,1:(size(X,2)-1),.1,[.1 1 10],k0);
        model.khat      = wRidge;
        model.fnlin     = nlfun;
        model.SDebars   = SDebars;
        model.rho       = rho;
end

if isfield(dspec.model, 'bilinearMode')
    covLabels = {dspec.covar(:).label};
    covDesc   = {dspec.covar(:).desc};
    switch upper(dspec.model.bilinearMode)
        case 'ON'
            brank = dspec.model.bilinearRank;
            
            % identify which coefficients are related to stimulus kernel
            [~,id] = findFile(covDesc, dspec.model.bilinearCovariate);
            grpCovLabels = covLabels(id);
            dspec.model.bilinearGroupList = id;
            nbil = length(dspec.model.bilinearGroupList); % number of coherence filters
            nprsperf = [dspec.covar(id).edim];
            assert(all(nprsperf==nprsperf(1)), 'all filters in bilinear group must be the same size')
            nprsperf = nprsperf(end); % # of params per coh filter
            ncohprs = nbil*nprsperf;
            iicoh = cell2mat(getDesignMatrixColIndices(dspec, grpCovLabels));
            kcohdims = [nprsperf,nbil];  % size of coherence-related filter kcoh
            
            % Do bilinear optimization (coordinate ascent)
            [khat, SDebars] = bilinearMixRegress_Poisson(X(ndx,:), Y(ndx), kcohdims, brank, ...
                iicoh, nlfun, dt);
            
            model.khat      = khat;
            model.fnlin     = nlfun;
            model.SDebars   = SDebars;
            
            return
        case 'OFF'
        otherwise
            error('Unknown bilinear mode');
    end
end

%% get Design matrix column indices
function [idx] = getDesignMatrixColIndices(dspec, covarLabels)
% Input
%   dpsec: design specification structure
%   covarLabels: 'str' or {'str'} - label(s) of the covariates
% Outut
%   idx: {} - column indices of the design matrix that correspond to the
%	    specified covariates

subIdxs = getGroupIndicesFromDesignSpec(dspec);

if ~iscell(covarLabels)
    covarLabels = {covarLabels};
end

labels={dspec.covar.label}';
labels=[labels num2cell((1:numel(labels))')]';
idxmap=struct(labels{:});
idx = cell(numel(covarLabels), 1);

for k = 1:numel(covarLabels)
    idx{k} = subIdxs{idxmap.(covarLabels{k})}(:);
end



%% get group indices
function subIdxs = getGroupIndicesFromDesignSpec(dspec)
% Cell of column indices that corresponds to each covariate in the design matrix
% subIdxs = getGroupIndicesFromDesignSpec(dspec)

subIdxs = {};
k = 0;

for kCov = 1:numel(dspec.covar)
    edim = dspec.covar(kCov).edim;
    subIdxs{kCov} = k + (1:edim); %#ok<AGROW>
    k = k + edim;
end










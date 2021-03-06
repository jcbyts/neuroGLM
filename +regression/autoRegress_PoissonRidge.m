function [wRidge,rho,SDebars,postHess,logevid] = autoRegress_PoissonRidge(xx,yy,fnlin,iirdge,rhoNull,rhovals,w0, dtbin)
% [wRidge,rho,SDebars,postHess,logevid] = autoRegress_PoissonRidge(xx,yy,fnlin,iirdge,rhoNull,rhovals,w0, dtbin)
%
% Computes empirical Bayes ridge regression filter estimate under a Poisson-GLM.
%
% Inputs:
%        xx [n x p] - stimulus (regressors)
%        yy [n x 1] - spike counts (response)
%     fnlin [1 x 1] - func handle for nonlinearity (must return func + 1st & 2nd derivs)
%     iirdg [q x 1] - weight indices for which to learn ridge parameter (columns of xx)
%   rhoNull [1 x 1] - fixed prior precision for other columns of xx
%   rhovals [v x 1] - list of prior precisions to use for grid search
%        w0 [p x 1] - initial estimate of regression weights (optional)
%
% Outputs:
%   wRidge [p x 1] - empirical bayes estimate of weight vector
%      rho [1 x 1] - estimate for prior precision (inverse variance)
%  SDebars [p x 1] - 1 SD error bars from posterior
% postHess [p x p] - Hessian (2nd derivs) of negative log-posterior at
%                    maximum ( inverse of posterior covariance)
%  logevid [1 x 1] - log-evidence for hyperparameters
%
% $Id$

import regression.*

% initialize filter estimate with MAP regression estimate, if necessary
if nargin == 6
    nw = size(xx,2); % length of stimulus filter
    rho0 = 5;        % initial value of ridge parameter
    Lprior = speye(nw)*rho0;
    Lprior(iirdge,iirdge) = rho0;
    XX = xx'*xx;
    XY = xx'*yy;
    w0 = (XX + Lprior) \ XY;
end

if ~exist('dtbin', 'var')
    dtbin=1;
end

% --- set prior and log-likelihood function pointers ---
mstruct.neglogli = @neglogli_poiss;
mstruct.logprior = @logprior_ridge;
mstruct.liargs = {xx,yy,fnlin,dtbin};
mstruct.priargs = {iirdge,rhoNull};

% --- Do grid search over ridge parameter -----
[hprsMax,wmapMax] = regression.gridsearch_GLMevidence(w0,mstruct,rhovals);
fprintf('best grid point: rho (precision)=%.1f\n', hprsMax);

% --- Do gradient ascent on evidence ----
[wRidge,rho,logevid,postHess] = findEBestimate_GLM(wmapMax,hprsMax,mstruct);

if nargout >2
    SDebars = sqrt(diag(inv(postHess)));
end

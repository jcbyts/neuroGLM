%% test bilinear poisson regression


ntb  = 15;
nsb  = 10;
nlin = 20;

wspace = exp(- ((1:nsb)-nsb/2).^2/15);
wtime  = cos( .5*(1:ntb));
wlin = exp(-cos(1:nlin));
wst = wtime(:)*wspace(:)';

figure(1); clf
subplot(2,2,1)
plot(wlin)
subplot(2,2,2)
imagesc(wst)


wtrue = [wlin(:); wst(:)]/10;
nk       = numel(wtrue); 
nsamples = 500;
X = randn(nsamples, nk);

y = poissrnd(exp(X*wtrue));

%%
wDims = [ntb nsb];

b=glmfit(X,y, 'poisson');
subplot(2,2,1)
imagesc(reshape(b(1:prod(wDims)), wDims))
title('glmfit')

%% linear bilinear regression

xx = X'*X;
xy = X'*y;

p = 1; % rank 1
indsbilin = (nlin+1):nk;
[wMLE,~,wt,wx,~] = regression.bilinearMixRegress_Poisson(X,y,wDims,p,indsbilin);
subplot(2,2,3)
imagesc(wt(:)*wx(:)')
title('bilinear MLE')

%% test with autoregress poisson
mstruct.nlfun=@expfun;
mstruct.dtbin=1;
iirdge=1:nk;
rhoNull=.1;
rhovals=[.1 1 10];
[wRidge,rho,SDebars,postHess,logevid] = autoRegress_PoissonRidge(X,y,mstruct.nlfun,iirdge,rhoNull,rhovals,w_hat, mstruct.dtbin);

mstruct.neglogli=@neglogli_poiss;
mstruct.logprior=@logprior_ridge;
mstruct.hyperprs=rho;

[wRho,~,wt,wx,wlin] = bilinearMixRegress_Poisson(X,y,wDims,p,indsbilin,mstruct);
subplot(2,2,4)
imagesc(wt(:)*wx(:)')

%%
if sum(wt)<0
    wt=-wt;
    wx=-wx;
end

nfun=@(x) x./norm(x);
% figure(1); clf
subplot(1,3,1)
plot(nfun(wt)); hold on
plot(nfun(wtime), 'k')
subplot(1,3,2)
plot(nfun(wx)); hold on
plot(nfun(wspace), 'k')
subplot(1,3,3)
plot(nfun(wtrue), 'k'); hold on
plot(nfun(w_hat))
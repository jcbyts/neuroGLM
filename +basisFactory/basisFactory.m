function bs = basisFactory(basisName, binfun, duration)
% make basis functions
% bs = basisFactory(basisName, binfun, duration)
% Possibe Names
% 
% nonlinearly scaled raised-cosine basis functions extra ones for refractory 
% bs = makeSpikeBasis(nBasis, binSize, rnge, stretch, nDelta)
% Possible inputs:
% nBasis, 
% binSize, 
% rnge, 
% stretch, 
% nDelta

switch basisName
    case {'history', 'spike'}
        if ~exist('duration', 'var'), duration=100;end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 8, binfun, 2);
        nDelta=1;
        bs.tr = [bs.tr(:,1:nDelta) bs.tr];
        bs.B = [[eye(nDelta); zeros(size(bs.B,1)-nDelta,nDelta)] bs.B];
        bs.edim = bs.edim + nDelta;
        
    case 'MTLIPcoupling'
        if ~exist('duration', 'var'), duration=100;end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 8, binfun, 10);
        
    case 'MTcoupling'
        if ~exist('duration', 'var'), duration=100;end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 8, binfun, 2);

    case 'LIPcoupling'
        if ~exist('duration', 'var'), duration=100;end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 8, binfun, 2);
        
    case 'LIPpulse'
        if ~exist('duration', 'var'), duration=[50 1e3];end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 10, binfun, 400);
        
    case 'MTpulse'
        if ~exist('duration', 'var'), duration=[50 400];end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 10, binfun, 250);
        
    case 'MTmotion'
        if ~exist('duration', 'var'), duration=[50 1000];end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 15, binfun, 500);
        
    case 'LIPmotion'
        if ~exist('duration', 'var'), duration=[50 1200];end
        bs=basisFactory.makeSmoothTemporalBasis('nonlinearly scaled cosine', duration, 15, binfun, 500);
        
    otherwise
        error('basis name not recognized')
        
end

bs.normalize;
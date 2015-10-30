function bases = makeSmoothTemporalBasis(shape, duration, nBases, binfun, nlOffset)
% Make Smooth Temporal Basis
% bases = makeSmoothTemporalBasis(shape, duration, nBases, binfun, nlOffset)
%
% Input
%   shape: 'raised cosine' or 'boxcar' or 'nonlinearly scaled cosine'
%   duration: the time that needs to be covered
%   nBases: number of basis vectors to fill the duration
%   binfun:
%
% Output
%   BBstm: basis object

if nargin < 1
    help basisFactory.makeSmoothTemporalBasis
    bases=[];
    return
end

nkbins = binfun(duration); % number of bins for the basis functions
ttb = repmat((1:nkbins)', 1, nBases); % time indices for basis

switch shape
    case 'nonlinearly scaled cosine'
        
        nlin = @(x)(log(x + 1e-20));
        invnl = @(x)(exp(x) - 1e-20);
        
        if nlOffset < 0
            error('nlOffset must be greater than 0');
        end
        if numel(duration)==2
            yrnge = nlin(duration + nlOffset);
        else
            yrnge = nlin([0 duration] + nlOffset);
        end
        binSize=1e3/binfun(1e3);
        db = diff(yrnge) / (nBases-1); % spacing between raised cosine peaks
        ctrs = yrnge(1):db:yrnge(2); % centers for basis vectors
        mxt = invnl(yrnge(2)+2*db) - nlOffset; % maximum time bin
        iht = (0:binSize:mxt)';
        ttb=repmat(iht/binSize, 1, nBases)+1;
        ff = @(x,c,dc) (cos(max(-pi, min(pi, (x-c)*pi/dc/2))) + 1)/2;
        BBstm = ff(repmat(nlin(iht + nlOffset), 1, nBases), repmat(ctrs, numel(iht), 1), db);
        bcenters = invnl(ctrs);
        
    case 'raised cosine'
        %   ^
        %  / \
        % /   \______
        %      ^
        %     / \
        % ___/   \___
        %         ^
        %        / \
        % ______/   \
        % For raised cosine, the spacing between the centers must be 1/4 of the
        % width of the cosine        
        dbcenter = nkbins / (3 + nBases); % spacing between bumps
        width = 4 * dbcenter; % width of each bump
        bcenters = 2 * dbcenter + dbcenter*(0:nBases-1);
        % location of each bump centers
        % bcenters = round(bcenters); % round to closest ms <-- BAD!!! NEVER DO THIS
        bfun = @(x,period)((abs(x/period)<0.5).*(cos(x*2*pi/period)*.5+.5));
        BBstm = bfun(ttb-repmat(bcenters,nkbins,1), width);
    case 'boxcar'
        width = nkbins / nBases;
        BBstm = zeros(size(ttb));
        bcenters = width * (1:nBases) - width/2;
        for k = 1:nBases
            idx = ttb(:, k) > ceil(width * (k-1)) & ttb(:, k) <= ceil(width * k);
            BBstm(idx, k) = 1 / sum(idx);
        end
    otherwise
        error('Unknown basis shape');
end

bases=Basis();
bases.type = [shape '@' mfilename];
bases.shape = shape;
bases.duration = duration;
bases.nBases = nBases;
bases.binfun = binfun;
bases.B = BBstm;
bases.edim = size(bases.B, 2);
bases.tr = ttb - 1;
bases.centers = bcenters;

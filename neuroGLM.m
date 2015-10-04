classdef neuroGLM < handle
    
    properties
        covar
        idxmap
        dm
        edim
        binfun
        binSize
        param
    end
    
    methods
        function obj=neuroGLM(varargin)
            obj.covar = struct();
            obj.idxmap = struct(); % reverse positional indexing
            
            if nargin > 0
                obj.initExperiment(varargin{:})
            end
        end
        
        %% initialize Experiment
        function initExperiment(obj, unitOfTime, binSize, param)
            % Initialize the Experiment structure that holds raw data
            % expt = initExperiment(unitOfTime, binSize, param)
            %
            %   unitOfTime: 'string' - 's' or 'ms' indicating a global unit of time
            %   binSize: [1] - duration of each time bin in units of unitOfTime
            %   param: (any) optional - any potentially useful extra information
            %       associated with the entire experiment record (experimental
            %       parameters)
            if nargin <=1
                help initExperiment
                return
            end
            
            assert(binSize > 0);
            assert(ischar(unitOfTime), 'Put a string for unit');
            
            obj.binSize  = binSize;
            obj.binfun   = @(t) (t == 0) + ceil(t/obj.binSize);
            
            
            if nargin > 3
                obj.param = param;
            else
                obj.param = struct();
            end
            
            obj.param.timeunit = unitOfTime;
            
        end
        
        
        %% remove constant columns
        function removeConstantCols(obj)
            % Remove columns with constant values (redandunt with the bias term)
            % dm = removeConstantCols(dm);
            %
            % If you want to also z-score, do it after removing constant columns
            
            if isfield(obj.dm, 'constCols')
                warning('This design matrix has already been cleared');
            end
            
            
            if isfield(obj.dm, 'zscore')
                error('Please zscore after removing constant columns');
            end
            
            v = var(obj.dm.X);
            
            obj.dm.constCols = v == 0;
            obj.dm.X = obj.dm.X(:, ~obj.dm.constCols);
        end
        
        %% z score design matrix
        function zscoreDesignMatrix(obj, colIndices)
            % z-scores each column of the design matrix. Helps when ill-conditioned.
            % dm = zscoreDesignMatrix(dm, colIndices);
            %
            % z-scoring removes the mean, and divides by the standard deviation.
            % CAUTION: this results in a non-sparse matrix if not used carefully!
            %
            % Input
            %   dm: design matrix structure (see compileSparseDesignMatrix)
            %   colIndices: [n x 1] optional - indicies of columns you want to z-score
            %	use getDesignMatrixColIndices to get indices from covariate labels
            %
            % Output
            %   dm: design matrix with z-scored columns
            %	it has the meta data to correctly reconstruct weights later
            
            if isfield(obj.dm, 'zscore')
                warning('This design matrix is already z-scored! Skipping...');
                return
            end
            
            if nargin == 1
                [obj.dm.X, obj.dm.zscore.mu, obj.dm.zscore.sigma] = zscore(obj.dm.X);
            else
                [X, zmu, zsigma] = zscore(obj.dm.X(:, colIndices));
                obj.dm.X(:, colIndices) = X;
                
                obj.dm.zscore.mu = zeros(size(obj.dm.X, 2), 1); % note that mean is not really zero
                obj.dm.zscore.mu(colIndices) = zmu;
                obj.dm.zscore.sigma = ones(size(obj.dm.X, 2), 1); % likewise
                obj.dm.zscore.sigma(colIndices) = zsigma;
            end
            
        end
        %% add Covariate with specified requirements
        function addCovariate(obj, trial, covLabel, desc, stimHandle, basisObj, offset, cond, plotOpts)
            % add covariate helper function
            % addCovariate(obj, trial, covLabel, desc, stimHandle, basisStruct, offset, cond, plotOpts)
            %
            % see also: addCovariateRaw, addCovariateSpikeTrain, addCovariateBoxcar, addCovariateTiming
            if nargin <= 1
                help addCovariate
                return
            end
            
            if ~ischar(covLabel)
                error('Covariate label must be a string');
            end
            
            if nargin < 4; desc = covLabel; end
            if nargin < 5; stimHandle = []; end
            
            if isfield(obj.idxmap, covLabel)
                warning('Label already added as a covariate');
                return
            end
            
            newIdx = numel(fieldnames(obj.idxmap)) + 1;
            
            if newIdx==1
                obj.covar=Covar(covLabel, desc, stimHandle);
            else
                obj.covar(newIdx)=Covar(covLabel, desc, stimHandle);
            end
            
            obj.idxmap.(covLabel) = newIdx;
            
            sdim = size(stimHandle(trial(1)), 2);
            obj.covar(newIdx).sdim = sdim;
            
            if nargin >= 6
                if isa(basisObj, 'Basis')
                    obj.covar(newIdx).basis = basisObj;
                    obj.covar(newIdx).edim  = basisObj.edim * sdim;
                else
                    error('Basis structure should be a basisObject');
                end
            else
                obj.covar(newIdx).edim = sdim;
            end
            
            if nargin >= 7
                obj.covar(newIdx).offset = obj.binfun(offset);
            else
                obj.covar(newIdx).offset = 0;
            end
            
            if nargin >= 8
                if ~isempty(cond) && ~isa(cond, 'function_handle')
                    error('Condition must be a function handle that takes trial');
                end
                obj.covar(newIdx).cond = cond;
            else
                obj.covar(newIdx).cond = [];
            end
            
            if nargin >= 9
                obj.covar(newIdx).plotOpts = plotOpts;
            end
            
            obj.edim = sum([obj.covar(:).edim]);
        end
        
        %% add Boxcar covariate
        function addCovariateBoxcar(obj, trial, covLabel, startLabel, endLabel, desc, varargin)
            % add boxcar covariate to model
            % addCovariateBoxcar(obj, trial, covLabel, startLabel, endLabel, desc, varargin)
            %
            % varargin: offset, cond, plotOpts
            %
            % see also: addCovariateRaw, addCovariateSpikeTrain, addCovariate, addCovariateTiming
            if nargin < 5; desc = covLabel; end
            
            assert(ischar(desc), 'Description must be a string');
            
            stimHandle = @(trial) basisFactory.boxcarStim(obj.binfun(trial.(startLabel)), obj.binfun(trial.(endLabel)), obj.binfun(trial.duration));
            
            obj.addCovariate(trial, covLabel, desc, stimHandle, varargin{:});
        end
        
        %% add raw continuous covariate (no basis)
        function addCovariateRaw(obj, trial, covLabel, desc, varargin)
            % Add the continuous covariate without basis function (instantaneous rel)
            %
            % varargin: offset, cond, plotOpts
            %
            
            if nargin < 3; desc = covLabel; end
            
            % assert(ischar(desc), 'Description must be a string');
            
            obj.addCovariate(trial, covLabel, desc, basisFactory.rawStim(covLabel), varargin{:});
            
        end
        
        %% add Spike train covariate
        function addCovariateSpiketrain(obj, trial, covLabel, stimLabel, desc, bs, varargin)
            % add spike train as covariate
            %
            % varargin: offset, cond, plotOpts
            %
            % addCovariateSpiketrain(trial, covLabel, stimLabel, desc, basisStruct, varargin)
            
            if nargin<=1
                help addCovariateSpiketrain
                return
            end
            
            if nargin < 5; desc = covLabel; end
            
            if nargin < 6
                bs = basisFactory.basisFactory('history', obj.binfun);
            end
            
            assert(ischar(desc), 'Description must be a string');
            
            offset = obj.binSize; % Make sure to be causal. No instantaneous interaction allowed.
            
            if numel(varargin) > 0 && isa(varargin{1}, 'function_handle')
                stimHandle = varargin{1};
                options={};
                if numel(varargin)>1
                    options=varargin(2:end);
                end
            else
                stimHandle = @(trial) basisFactory.deltaStim(obj.binfun(trial.(stimLabel)), obj.binfun(trial.duration));
%                 stimHandle = @(trial) basisFactory.deltaStim(obj.binfun(trial.(stimLabel)+obj.binSize), obj.binfun(trial.duration));
                options = varargin;
            end
            
            obj.addCovariate(trial, covLabel, desc, stimHandle, bs, offset, options{:});
        end
        
        %% add Timing Covariate
        function addCovariateTiming(obj, trial, covLabel, stimLabel, desc, varargin)
            % Add a timing covariate based on the stimLabel.
            %
            % varargin: offset, cond, plotOpts
            %
            % addCovariateTiming(trial, covLabel, stimLabel, desc, varargin)
            
            if nargin < 4; stimLabel = covLabel; end
            if nargin < 5; desc = covLabel; end
            
            if isempty(stimLabel)
                stimLabel = covLabel;
            end
            
            if isempty(desc)
                desc=stimLabel;
            end
                        
            stimHandle = @(trial) basisFactory.deltaStim(obj.binfun(trial.(stimLabel)), obj.binfun(trial.duration(1)));
            
            obj.addCovariate(trial, covLabel, desc, stimHandle, varargin{:});
            
        end
        
        %% combine fitted weights
        function [wout] = combineWeights(obj, w)
            % Combine the weights per column in the design matrix per covariate
            % [wout] = combineWeights(weightVector)
            %
            % Input
            %   dm: design matrix structure
            %   w: weight on the basis functions
            %
            % Output
            %   wout.(label).data = combined weights
            %   wout.(label).tr = time axis
            
            obj.binSize;
            wout = struct();
            
            if isfield(obj.dm, 'biasCol') % undo z-score operation
                if isfield(obj.dm, 'zscore') && numel(obj.dm.zscore.mu) == numel(w) % remove bias from zscore
                    zmu  = obj.dm.zscore.sigma(obj.dm.biasCol);
                    zsig = obj.dm.zscore.mu(obj.dm.biasCol);
                    obj.dm.zscore.sigma(obj.dm.biasCol) = [];
                    obj.dm.zscore.mu(obj.dm.biasCol) = [];
                else
                    zmu  = 0;
                    zsig = 1;
                end
                wout.bias = w(obj.dm.biasCol)*zsig + zmu; % un-z-transform the bias
                w(obj.dm.biasCol) = [];
            end
            
            if isfield(obj.dm, 'zscore') % undo z-score operation
                w = (w .* obj.dm.zscore.sigma(:)) + obj.dm.zscore.mu(:);
            end
            
            if isfield(obj.dm, 'constCols') % put back the constant columns
                w2 = zeros(obj.edim, 1);
                w2(~obj.dm.constCols) = w; % first term is bias
                w = w2;
            end
            
            if numel(w) ~= obj.edim
                error('Expecting w to be %d dimension but it''s [%d]', ...
                    obj.edim, numel(w));
            end
            
            startIdx = [1 (cumsum([obj.covar(:).edim]) + 1)];
            
            for kCov = 1:numel(obj.covar)
                bs = obj.covar(kCov).basis;
                
                if isempty(bs)
                    w_sub = w(startIdx(kCov) + (1:obj.covar(kCov).edim) - 1);
                    wout.(obj.covar(kCov).label).tr = ((1:size(w_sub, 1))-1 + obj.covar(kCov).offset) * obj.binSize;
                    wout.(obj.covar(kCov).label).data = w_sub;
                    continue;
                end
                
                assert(isa(bs, 'Basis'), 'Basis is not a Basis?');
                
                sdim = obj.covar(kCov).edim / bs.edim;
                wout.(obj.covar(kCov).label).data = zeros(size(bs.B, 1), sdim);
                for sIdx = 1:sdim
                    w_sub = w(startIdx(kCov) + (1:bs.edim)-1 + bs.edim * (sIdx - 1));
                    w2_sub = sum(bsxfun(@times, bs.B, w_sub(:)'), 2);
                    wout.(obj.covar(kCov).label).data(:, sIdx) = w2_sub;
                end
                wout.(obj.covar(kCov).label).tr = ...
                    ((bs.tr(:, 1) + obj.covar(kCov).offset) * obj.binSize) * ones(1, sdim);
            end
            
        end
        
        %% Add Bias Column to design matrix
        function addBiasColumn(obj,flag)
            % Add a column of ones as the first (or last) column to estimate the bias
            % (DC term)
            % addBiasColumn(flag);
            %
            % Some regression packages do not allow separate bias estimation, and
            % a constant column of ones must be added to the design matrix. The weight
            % that corresponds to this column will be the bias later.
            %
            % To put the bias column on the right of the design matrix, call:
            % addBiasColumn('right')
            
            if isfield(obj.dm, 'biasCol')
                warning('This design matrix already has a bias column! Skipping...');
                return
            end
            
            if nargin < 2
                flag = 'left'; % bias column defaults to the left of the design matrix
            end
            
            switch flag
                case 'left'
                    obj.dm.X = [ones(size(obj.dm.X, 1), 1), obj.dm.X];
                    obj.dm.biasCol = 1; % indicating that the first column is the bias column
                case 'right'
                    n = obj.edim+1;
                    obj.dm.X = [obj.dm.X, ones(size(obj.dm.X, 1), 1)];
                    obj.dm.biasCol = n; % indicating that the last column is the bias column
            end
            
        end
        
        %% Compile Design Matrix
        function compileDesignMatrix(obj, trial, trialIndices)
            % Compile information from experiment according to given DesignSpec
            if ~isempty(obj.dm) && (numel(trialIndices)==numel(obj.dm.trialIndices)) && all(sum(bsxfun(@eq, trialIndices(:),obj.dm.trialIndices),2))
                return
            end
            
            fprintf('building the design matrix\n')
            obj.dm=struct();
            subIdxs = obj.getGroupIndicesFromDesignSpec;
            
            trialT = obj.binfun([trial(trialIndices).duration]);
            totalT = sum(trialT);
            X      = zeros(totalT, obj.edim);
            
            for k = 1:numel(trialIndices)
                kTrial = trialIndices(k);
                ndx = (sum(trialT(1:k))-(trialT(k)-1)):sum(trialT(1:k));
                
                for kCov = 1:numel(obj.covar) % for each covariate
                    
                    sidx = subIdxs{kCov};
                    
                    if ~isempty(obj.covar(kCov).cond) && ~obj.covar(kCov).cond(trial(kTrial))
                        continue;
                    end
                    
                    stim = obj.covar(kCov).stim(trial(kTrial)); % either dense or sparse
                    stim = full(stim);
                    
                    if ~isempty(obj.covar(kCov).basis)
                        X(ndx, sidx) = obj.covar(kCov).basis.convolve(stim, obj.covar(kCov).offset);
%                         X(ndx, sidx) = basisFactory.convBasis(stim, obj.covar(kCov).basis, obj.covar(kCov).offset);
                    else
                        X(ndx, sidx) = stim;
                    end
                end
            end
            
            obj.dm.X = X;
            obj.dm.trialIndices = trialIndices;
            
            % Check sanity of the design
            if any(~isfinite(obj.dm.X(:)))
                warning('Design matrix contains NaN or Inf...this is not good!');
            end
            
        end
        
        %% get binned spike train
        function y = getBinnedSpikeTrain(obj, trial, spLabel, trialIdx)
            % y: a sparse column vector representing the concatenated spike trains
%             getBinnedSpikeTrain(obj, trial, spLabel, trialIdx)
            
            sts = cell(numel(trialIdx), 1);
            endTrialIndices = [0 cumsum(obj.binfun([trial(trialIdx).duration]))];
            nT = endTrialIndices(end); % how many bins total?
            
            for k = 1:numel(trialIdx)
                kTrial = trialIdx(k);
                bst = endTrialIndices(k) + obj.binfun(trial(kTrial).(spLabel));
                sts{k} = bst(:);
            end
            
            sts = cell2mat(sts);
            sts(sts>nT) = [];
            y = sparse(sts, 1, 1, nT, 1);
            
        end
        
        %% get Design matrix column indices
        function [idx] = getDesignMatrixColIndices(obj, covarLabels)
            % Input
            %   dpsec: design specification structure
            %   covarLabels: 'str' or {'str'} - label(s) of the covariates
            % Outut
            %   idx: {} - column indices of the design matrix that correspond to the
            %	    specified covariates
            
            subIdxs = obj.getGroupIndicesFromDesignSpec();
            
            if ~iscell(covarLabels)
                covarLabels = {covarLabels};
            end
            
            idx = cell(numel(covarLabels), 1);
            
            for k = 1:numel(covarLabels)
                idx{k} = subIdxs{obj.idxmap.(covarLabels{k})}(:);
            end
            
        end
        
        %% get group indices
        function subIdxs = getGroupIndicesFromDesignSpec(obj)
            % Cell of column indices that corresponds to each covariate in the design matrix
            % subIdxs = getGroupIndicesFromDesignSpec(dspec)
            
            subIdxs = {};
            k = 0;
            
            for kCov = 1:numel(obj.covar)
                ddim = obj.covar(kCov).edim;
                subIdxs{kCov} = k + (1:ddim); %#ok<AGROW>
                k = k + ddim;
            end
            
        end
        
        %% get continuous response variable
        function y = getResponseVariable(obj, trial, label, trialIdx)
            % y: a column vector or matrix representing the concatenated continuous variable
            
            if ~isempty(trial)
                error('First argument must be an experiment structure');
            end
            
            if ~isfield(obj.expt.type, label)
                error('Label [%s] is not registered in the experiment structure', label);
            end
            
            % Check that the label corresponds to a continuous variable
            if ~strcmp(obj.expt.type.(label), 'continuous')
                error('Type of label [%s] is not continuous', label);
            end
            
            % put everything in a cell
            ycell = cell(numel(trialIdx), 1);
            
            for kTrial = trialIdx(:)'
                ycell{kTrial} = trial(kTrial).(label);
            end
            
            y = cell2mat(ycell);
        end
        
        %% get covariate timing
        function X = getCovariateTiming(obj, trial, trialIndices)
            % Compile information from experiment according to given DesignSpec
            % X = getCovariateTiming(trial,trialIndices)
            
            trialT = ceil([trial(trialIndices).duration]/obj.binSize);
            totalT = sum(trialT);
            
            nCovariates = numel(obj.covar);
            
            X      = zeros(totalT, nCovariates);
            
            for k = 1:numel(trialIndices)
                kTrial = trialIndices(k);
                ndx = sum(trialT(1:k))-(trialT(k)-1):sum(trialT(1:k));
                
                for kCov = 1:nCovariates
                    if ~isempty(obj.covar(kCov).cond) && ~obj.covar(kCov).cond(trial(kTrial))
                        continue;
                    end
                    
                    stim = obj.covar(kCov).stim(trial(kTrial)); % either dense or sparse
                    X(ndx, kCov) = full(stim);
                end
            end
        end
        
        
    end
end







% % %
% %  if nargin<2
% %                 fprintf(['*****************************************************************************\n', ...
% %                     'Help for addCovariate\n', ...
% %                     '*****************************************************************************\n', ...
% %                     'addCovariate(obj, covLabel, desc, stimHandle, basisStruct, offset, cond, plotOpts)\n', ...
% %                     'Add a covariate to the design specification object.\n'...
% %                     'obj.dspec = addCovariate(dspec, covLabel, desc, stimHandle, basisStruct, offset, cond, plotOpts)\n',...
% %                     'Input:\n'...
% %                     '  dspec: \tstruct - design specification object (see initDesignSpec)\n'...
% %                     '  covLabel: \tstring - name to refer to the covariate\n'...
% %                     '  desc: \tstring - human readable description of the covariate\n'...
% %                     '  stimHandle: @(trial, expt) -> [T x m] the raw stimulus design before\n'...
% %                     '  applying the basis functions.\n\n'...
% %                     '  basisStruct: struct - tempora,l basis functions that will be convolved with,\n'...
% %                     '  the output of the stimHandle. See +basisFactory functions \n(e.g.'...
% %                     '  +basisFactory.makeSmoothTempo,ralBasis),\n'...
% %                     '  offset: \t[1] optional/default:,0 - number of **time bins** to shift the '...
% %                     '  regressors.\n Negative (positive) integers represent acausal (causal) effects.\n'...
% %                     '  cond: @(,/nt...rial) -> boolean optional: condition for which the covariate will,\n'...
% %                     '  be included. For example, if only trials where ''choice'' is 1 to include\n',...
% %                     '  the current covariate, use @(trial) (trial.choice == 1),\n', ...
% %                     'Output:\n',...
% %                     '   dspec: updated design specification object\n',...
% %                     'See also: addCovariateTiming, addCovariateRaw, addCovariateBoxcar, addCovariateSpiketrain\n'])
% %             end
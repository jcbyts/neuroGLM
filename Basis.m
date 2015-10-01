classdef Basis < handle
    properties
        type
        shape
        duration
        nBases
        binfun
        B
        edim
        tr
        centers
        normalized
        orthogonalized
    end
    
    methods
        function b=Basis()
            b.type='none';
            b.shape='none';
            b.duration=1;
            b.nBases=1;
            b.B=1;
            b.tr=1;
            b.edim=1;
            b.centers=0;
            b.normalized=false;
            b.orthogonalized=false;
        end
        
        function plot(b)
            plot(b.tr, b.B)
        end
        
        function normalize(b)
            b.B = bsxfun(@rdivide, b.B, sum(b.B));
            b.normalized=true;
        end

        function orthogonalize(b)
            b.B = orth(b.B);
            b.orthogonalized=true;
        end

        function X = convolve(b, stim, offset)
            % Convolve basis functions to the covariate matrix
            %
            %   b.convolve(stim, offset)
            %
            %     stim: [T x dx]     - covariates over time
            %   offset: [1] optional - shift in time
            
            if nargin < 3
                offset = 0;
            end
            
            [~, dx] = size(stim);
            
            % zero pad stim to account for the offset
            if offset < 0 % anti-causal
                stim = [stim; zeros(-offset, dx)];
            elseif offset > 0; % push to future
                stim = [zeros(offset, dx); stim];
            end
            
            if issparse(stim) || nnz(stim) < 20;
                X = basisFactory.temporalBases_sparse(stim, b.B);
            else
                X = basisFactory.temporalBases_dense(stim, b.B);
            end
            
            if offset < 0 % anti-causal
                X = X(-offset+1:end, :);
            elseif offset > 0
                X = X(1:end-offset, :);
            end
        end
        
    end
end
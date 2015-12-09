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
        
        function S=saveobj(b)
            c=b;
%             c=funh2struct(b, false);
            S=struct();
            
            S.type=c.type;
            S.shape=c.shape;
            S.duration=c.duration;
            S.nBases=c.nBases;
            S.binfun=c.binfun;
            S.B=c.B;
            S.edim=c.edim;
            S.tr=c.tr;
            S.centers=c.centers;
            S.normalized=c.normalized;
            S.orthogonalized=c.orthogonalized;
        end
        
        
    end
    
    methods(Static)
        
        function obj=loadobj(s)
            if isstruct(s)
                newObj=Basis;
                fields=fieldnames(s);
                for kField=1:numel(fields)
                    newObj.(fields{kField})=s.(fields{kField});
                end
            else
                newObj=s;
            end
%             obj=struct2funh(newObj,false);
            obj=newObj;
        end
    end
end
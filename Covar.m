classdef Covar < handle
    
    properties
        label
        desc
        stim    % function_handle
        offset  % offset (in bins)
        cond    % function_handle (condition that must be met to add)
        basis   % basis
        edim
        sdim
    end
    
    methods
        
        function c=Covar(label, desc, stim)  % constructor
            % Covariate Class
            % c=Covar(label, desc, stim, offset, cond, basis)
            assert(isa(stim, 'function_handle'), 'covariate needs a function handle to get stimulus')
            assert(isa(label, 'char'), 'label must be a string')
            assert(isa(desc,  'char'), 'desc must be a string')
            
            c.label=label;
            c.desc=desc;
            c.stim=stim;
           
        end
        
        function plotBasis(c)
            c.basis.plot
            title(c.desc)
        end
%         
        function S=saveobj(c)
            d=c;
            if isa(d.stim, 'function_handle')
            f=functions(d.stim);
            if numel(f.workspace)>1
                if isfield(f.workspace{2}, 'obj')
                    f.workspace{2}=rmfield(f.workspace{2}, 'obj');
                end
                
                if isfield(f.workspace{2}, 'trial')
                    f.workspace{2}=rmfield(f.workspace{2}, 'trial');
                end
                
                if isfield(f.workspace{2}, 'stimHandle')
                    f.workspace{2}=rmfield(f.workspace{2}, 'stimHandle');
                end
                
                if isfield(f.workspace{2}, 'bf')
                    f.workspace{2}=rmfield(f.workspace{2}, 'bf');
                end
                
                if isfield(f.workspace{2}, 'varargin')
                    f.workspace{2}=rmfield(f.workspace{2}, 'varargin');
                end
            end
% %             d=funh2struct(c, false);
% %             struct2funh(c,false);
%             
            S.label=d.label;
            S.desc=d.desc;
%             if numel(d.stim.workspace)>1 && isfield(d.stim.workspace{2}, 'bs')
% %                 d.stim.workspace{1}.bs=d.stim.workspace{2}.bs;
% %                 if isfield(d.stim.workspace{1}, 'obj')
% %                     d.stim.workspace{1}=rmfield(d.stim.workspace{1}, 'obj');
% %                 end
% % %                 if isfield(d.stim.workspace{2}, 'binningMatrix')
% % %                     d.stim.workspace{1}.binningMatrix=d.stim.workspace{2}.binningMatrix;
% % %                 end
% %                 d.stim.workspace(2)=[];
% %             end
            S.stim=f;
            S.offset=d.offset;
            S.cond=d.cond;
            S.basis=d.basis;
            S.edim=d.edim;
            S.sdim=d.sdim;
        end
    end
    
    methods(Static)
        
        function obj=loadobj(s)
            d=struct2funh(s,false);
            obj=Covar(d.label, d.desc, d.stim);
            obj.offset=d.offset;
            obj.cond=d.cond;
            obj.basis=d.basis;
            obj.edim=d.edim;
            obj.sdim=d.sdim;
        end
    end
    
end
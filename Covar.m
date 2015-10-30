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
    end
    
    methods(Static)
        function sobj=saveobj(c)
            c=functionHandleRepair(c, false);
            sobj=struct();
            props=properties(c);
            for kProp=1:numel(properties)
                if isa(c.(props{kProp}), 'handle')
                    sobj.(props{kProp})=c.(props{kProp}).saveobj;
                else
                    sobj.(props{kProp})=c.(props{kProp});
                end
            end
            
        end
        
        function obj=loadobj(s)
            if isstruct(s)
                newObj=glmspike;
                fields=fieldnames(s);
                for kField=1:numel(fields)
                    newObj.(fields{kField})=s.(fields{k});
                end
            else
                newObj=s;
            end
            obj=functionHandleRepair(newObj,false);
        end
    end
    
end
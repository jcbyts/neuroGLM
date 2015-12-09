classdef supersaver < handle
    % handle class with modified save and load functions
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

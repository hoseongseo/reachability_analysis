classdef Link2D < DynSys
    properties
        xD % desired state
        uMax % Input
        dMax % Disturbance
        dims
    end
    
    methods
        function obj = Link2D(x, uMax, dMax, dims)
            
            if numel(x) ~= obj.nx
                error('Initial state does not have right dimension!');
            end
            
            if ~iscolumn(x)
                x = x';
            end
            
            % Basic vehicle properties
            obj.pdim = [find(dims == 1) find(dims == 2)]; % Position dimensions
            obj.hdim = find(dims == 3);
            obj.nx = length(dims);
            obj.nu = 3;
            obj.nd = 3;
            
            obj.x = x;
            obj.xhist = obj.x;
            
            obj.xD = x;
            obj.uMax = uMax;
            obj.dMax = dMax;
            obj.dims = dims;
        end
        
    end % end methods
end % end classdef

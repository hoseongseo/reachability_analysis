classdef HssQuad1 < DynSys
    properties
        xD % desired state
        uMax % Input
        dMax % Disturbance
        dims
        
        mass
        gains
        cD
    end
    
    methods
        function obj = HssQuad1(x, uMax, dMax, dims, mass, gains, cD)
            %%%%%
            % x(1) = x
            % x(2) = dx
            % x(3) = theta
            % u(1) = theta_d
            % mass : mass
            % gains = ktheta_p : inner controller gain
            % cD : drag coefficient
            %%%%%
            % dx(1) = x(2)
            % dx(2) = -cD/m * x(2) + 9.8*sin(x(3)) + D(1)
            % dx(3) = -ktheta_p * x(3) + ktheta_p * u(1) + D(2)
            
            if numel(x) ~= obj.nx
                error('Initial state does not have right dimension!');
            end
            
            if ~iscolumn(x)
                x = x';
            end
            
            % Basic vehicle properties
            obj.nx = length(dims);
            obj.nu = 1;
            obj.nd = 2;
            
            obj.x = x;
            obj.xhist = obj.x;
            
            obj.xD = x;
            obj.uMax = uMax;
            obj.dMax = dMax;
            obj.dims = dims;
            obj.mass = mass;
            obj.gains = gains;
            obj.cD = cD;
        end
        
    end % end methods
end % end classdef

classdef Subsys1 < DynSys
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
        function obj = Subsys1(x, uMax, dMax, dims, mass, gains, cD)
            %%%%%
            % x(1) = theta
            % x(2) = force
            % x(3) = x
            % x(4) = dx
            % u(1) = theta_d
            % u(2) = force_d
            % mass : mass
            % gains = [ktheta_p, kf_p] : inner controller gain
            % cD : drag coefficient
            %%%%%
            % dx(1) = -ktheta_p * x(1) + ktheta_p * u(1) + D(1)
            % dx(2) = -kf_p * x(2) + kf_p * u(2) + D(2)
            % dx(3) = x(4)
            % dx(4) = -cD/m * x(4) + 9.8*sin(x(1)) + x(2)/m*sin(x(1)) + D(3)
            
            if numel(x) ~= obj.nx
                error('Initial state does not have right dimension!');
            end
            
            if ~iscolumn(x)
                x = x';
            end
            
            % Basic vehicle properties
            obj.nx = length(dims);
            obj.nu = 2;
            obj.nd = 3;
            
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

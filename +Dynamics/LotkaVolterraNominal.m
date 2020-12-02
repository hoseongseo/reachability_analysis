classdef LotkaVolterraNominal < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
    end
    methods
        function this = LotkaVolterraNominal
            this = this@Dynamics.Base(2,1,1);
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            f_ = [3*x(1,:)*(1-x(2,:)); 3*x(2,:)*(x(1,:)-1)];
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            dfdx_ = [3*(1-x(2)), -3*x(1);...
                    3*x(2), 3*(x(1)-1)];
        end
    end
end
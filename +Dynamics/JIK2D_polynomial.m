classdef JIK2D_polynomial < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
    end
    properties (GetAccess = public, SetAccess = private)
    end
    methods
        function this = JIK2D_polynomial()
            this = this@Dynamics.Base(2,1,1);
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            f_ = [x(2,:); -(0.9930*x(1,:) - 0.1948*x(1,:).^3) - x(2,:)] + 0.102*[0; 1]*w;
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            dfdx_ = [0, 1; -cos(x(1)), -1];
        end
        function dfdw_ = dfdw(this, x, u, w, t)
            dfdw_ = 0.102*[0; 1];
        end
    end
end
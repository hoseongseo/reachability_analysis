classdef LotkaVolterra < SystemBase
    properties (GetAccess = private, SetAccess = private)
        wMax;
    end
    methods
        function this = LotkaVolterra(wMax)
            this = this@SystemBase(2,1,1);
            this.wMax = wMax;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            f_ = [(this.wMax*w+3)*x(1,:).*(1-x(2,:));...
                  (this.wMax*w+3)*x(2,:).*(x(1,:)-1)];
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            dfdx_ = [(this.wMax*w+3)*(1-x(2)), -(this.wMax*w+3)*x(1);...
                (this.wMax*w+3)*x(2), (this.wMax*w+3)*(x(1)-1)];
        end
        function dfdw_ = dfdw(this, x, u, w, t)
            dfdw_ = this.wMax*[x(1)*(1-x(2)); x(2)*(x(1)-1)];
        end
    end
end
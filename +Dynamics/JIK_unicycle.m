classdef JIK_unicycle < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
    end
    properties (GetAccess = public, SetAccess = private)
        vNominal;
        wNominal;
        wMax;
    end
    methods
        function this = JIK_unicycle(vNominal, wNominal, wMax)
            this = this@Dynamics.Base(3,1,1);
            this.vNominal = vNominal;
            this.wNominal = wNominal;
            this.wMax = wMax;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            f_ = [this.wNominal*x(2,:) + this.vNominal*(x(3,:).*x(3,:)*0.5);...
                -this.wNominal*x(1,:) + this.vNominal*(x(3,:) - (1/6)*x(3,:).*x(3,:).*x(3,:));...
                this.wMax*w*ones(1,size(x,2))];
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            dfdx_ = [0, this.wNominal, this.vNominal*sin(x(3));...
                -this.wNominal, 0, this.vNominal*cos(x(3));...
                0, 0, 0];
        end
        function dfdw_ = dfdw(this, x, u, w, t)
            dfdw_ = this.wMax*[0; 0; 1];
        end
    end
end
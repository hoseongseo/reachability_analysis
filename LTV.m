classdef LTV < SystemBase
    properties (GetAccess = private, SetAccess = private)
        A
        D
    end
    methods
        function this = LTV(nx, nu, nw, Amat, Dmat)
            this = this@SystemBase(nx,nu,nw);
            this.A = Amat;
            this.D = Dmat;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            f_ = this.A(t) * x + this.D(t) * w;
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            dfdx_ = this.A(t);
        end
        function dfdw_ = dfdw(this, x, u, w, t)
            dfdw_ = this.D(t);
        end
    end
end
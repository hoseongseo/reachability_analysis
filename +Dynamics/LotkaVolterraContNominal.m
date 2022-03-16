classdef LotkaVolterraContNominal < Dynamics.Base
    properties (GetAccess = public, SetAccess = private)
        xN
        tN
    end
    methods
        function this = LotkaVolterraContNominal(x0, t)
            this = this@Dynamics.Base(2,1,1);
            this.xN = this.rk4_forward(x0, zeros(this.Nu, length(t)), t);
            this.tN = t;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            f_ = [3*x(1,:)*((1-t)-x(2,:)); 3*x(2,:)*(x(1,:)-(1-t))];
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            dfdx_ = [3*((1-t)-x(2)), -3*x(1);...
                    3*x(2), 3*(x(1)-(1-t))];
        end
    end
    methods
        function xt = interpolate(this, tq)
            k1 = find((this.tN - tq) < 0, 1, 'last');
            k2 = find((this.tN - tq) > 0, 1, 'first');
            if (k2 - k1) == 1
                % Choice 1: linear interpolation
                % xt = this.xN(:,k1) + (this.xN(:,k2)-this.xN(:,k1))*(tq-this.tN(k1))/(this.tN(k2)-this.tN(k1));

                % Choice 2: local integration
                x0 = this.xN(:,k1);
                x_integ = this.rk4_forward(x0, zeros(this.Nu,2), [this.tN(k1), tq]);
                xt = x_integ(:,end);
                
                % Choice 3: piecewise constant assumption
                % xt = this.xN(:,k1);
            elseif (k2 - k1) == 2
                xt = this.xN(:,0.5*(k1+k2));
            elseif isempty(k1)
                xt = this.xN(:,1);
            elseif isempty(k2)
                xt = this.xN(:,end);
            end
        end
    end
end
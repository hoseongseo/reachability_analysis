classdef MultirotorNominal < Dynamics.Base
    properties (GetAccess = public, SetAccess = private)
        xN
        tN
    end
    methods
        function this = MultirotorNominal(x0, u, t)
            this = this@Dynamics.Base(8,3,3);
            this.xN = this.rk4_forward(x0, u, t);
            this.tN = t;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
%             f_ = [x(4,:);...
%                   x(5,:);...
%                   x(6,:);...
%                   u(1)*x(7,:) + w(1);...
%                   -u(1)*x(8,:) + w(2);...
%                   u(1) - 9.8 + w(3);...
%                   u(2)*ones(1,size(x,2));...
%                   u(3)*ones(1,size(x,2))];
            
            cp = 1 - 0.5*x(7,:).*x(7,:);
%             cp = 1;
%             sp = x(7,:) - (1/6)*(x(7,:).*x(7,:).*x(7,:));
            sp = x(7,:);
            ct = 1 - 0.5*x(8,:).*x(8,:);
%             ct = 1;
%             st = x(8,:) - (1/6)*(x(8,:).*x(8,:).*x(8,:));
            st = x(8,:);
            
            f_ = [x(4,:);...
                  x(5,:);...
                  x(6,:);...
                  u(1)*cp.*st + w(1);...
                  -u(1)*sp + w(2);...
                  u(1)*cp.*ct - 9.8 + w(3);...
                  u(2)*ones(1,size(x,2));...
                  u(3)*ones(1,size(x,2))];
        end
%         function dfdx_ = dfdx(this, x, u, w, t)
%             dfdx_ = zeros(this.Nx,this.Nx);
%             dfdx_(1,4) = 1.0;
%             dfdx_(2,5) = 1.0;
%             dfdx_(3,6) = 1.0;
%             dfdx_(4,7) = u(1);
%             dfdx_(5,8) = -u(1);
%         end
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
classdef LotkaVolterra < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
        wMax; % maximum disturbance
        sys0; % nominal dynamics
    end
    properties (GetAccess = public, SetAccess = private)
        xN;   % nominal trajectory
        tN;   % time of xN
    end
    methods
        function this = LotkaVolterra(wMax, sys0, x0, tN)
            this = this@Dynamics.Base(sys0.Nx,sys0.Nu,sys0.Nw);
            this.wMax = wMax;
            this.sys0 = sys0;
            this.xN = sys0.rk4_forward(x0, zeros(sys0.Nu,length(tN)), tN);
            this.tN = tN;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            q_ = this.interpolate( t );
            
            fq_ = this.sys0.f(q_, zeros(this.Nu,1), zeros(this.Nw,1), t);
            
            f_ = [(this.wMax*w+3)*(x(1,:)+q_(1)).*(1-(x(2,:)+q_(2)));...
                  (this.wMax*w+3)*(x(2,:)+q_(2)).*(x(1,:)+q_(1)-1)] - fq_;
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            q_ = this.interpolate( t );
            
            dfdx_ = [(this.wMax*w+3)*(1-(x(2)+q_(2))), -(this.wMax*w+3)*(x(1)+q_(1));...
                (this.wMax*w+3)*(x(2)+q_(2)), (this.wMax*w+3)*(x(1)+q_(1)-1)];
        end
        function dfdw_ = dfdw(this, x, u, w, t)
            q_ = this.interpolate( t );
            
            dfdw_ = this.wMax*[(x(1)+q_(1))*(1-(x(2)+q_(2))); (x(2)+q_(2))*(x(1)+q_(1)-1)];
        end
    end
    methods (Hidden)
        function xt = interpolate(this, tq)
            k1 = find((this.tN - tq) < 0, 1, 'last');
            k2 = find((this.tN - tq) > 0, 1, 'first');
            if (k2 - k1) == 1
                % Choice 1: linear interpolation
%                 xt = this.xN(:,k1) + (this.xN(:,k2)-this.xN(:,k1))*(tq-this.tN(k1))/(this.tN(k2)-this.tN(k1));

                % Choice 2: local integration
%                 x0 = this.xN(:,k1);
%                 xinteg = this.sys0.rk4_forward(x0, zeros(this.Nu,2), [this.tN(k1), tq]);
%                 xt = xinteg(:,end);
                
                % Choice 3: piecewise constant assump.
                xt = this.xN(:,k1);
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
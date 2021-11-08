classdef LotkaVolterra2 < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
        wMax; % maximum disturbance
    end
    properties (GetAccess = public, SetAccess = private)
%         xN;   % nominal trajectory
%         tN;   % time of xN
    end
    methods
        function this = LotkaVolterra2(wMax)
            this = this@Dynamics.Base(2,1,1);
            this.wMax = wMax;
%             this.sys0 = sys0;
%             this.xN = sys0.xN;
%             this.tN = sys0.tN;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
%             q_ = this.sys0.interpolate( t );
%             fq_ = this.sys0.f(q_, zeros(this.Nu,1), zeros(this.Nw,1), t);
            
            f_ = [(this.wMax*w+3)*x(1,:).*(1-x(2,:));...
                  (this.wMax*w+3)*x(2,:).*(x(1,:)-1)];
        end
        function dfdx_ = dfdx(this, x, u, w, t)
%             q_ = this.sys0.interpolate( t );
            
            dfdx_ = [(this.wMax*w+3)*(1-x(2)), -(this.wMax*w+3)*x(1);...
                (this.wMax*w+3)*x(2), (this.wMax*w+3)*(x(1)-1)];
        end
        function dfdw_ = dfdw(this, x, u, w, t)
%             q_ = this.sys0.interpolate( t );
            
            dfdw_ = this.wMax*[x(1)*(1-x(2)); x(2)*(x(1)-1)];
        end
    end
end
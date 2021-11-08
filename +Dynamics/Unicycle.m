classdef Unicycle < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
        wMax; % maximum disturbance
        sys0; % nominal dynamics
    end
    properties (GetAccess = public, SetAccess = private)
        xN;   % nominal trajectory
        tN;   % time of xN
    end
    methods
        function this = Unicycle(sys0, wMax)
            this = this@Dynamics.Base(sys0.Nx,sys0.Nu,sys0.Nw);
            this.wMax = wMax;
            this.sys0 = sys0;
            this.xN = sys0.xN;
            this.tN = sys0.tN;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            q_ = this.sys0.interpolate( t );
            
            fq_ = this.sys0.f(q_, zeros(this.Nu,1), zeros(this.Nw,1), t);
            
%             f_ = [ (1-x(3,:).^2*0.5)*cos(q_(3)) - (x(3,:) - (x(3,:).^3*(1/6)))*sin(q_(3));...
%                    (x(3,:) - (x(3,:).^3*(1/6)))*cos(q_(3)) + (1-x(3,:).^2*0.5)*sin(q_(3));...
%                    u + this.wMax*w*ones(1,size(x,2)) ] - fq_;
               
            f_ = [ (1-x(3,:).^2*0.5)*cos(q_(3)) - (x(3,:))*sin(q_(3));...
                   (x(3,:))*cos(q_(3)) + (1-x(3,:).^2*0.5)*sin(q_(3));...
                   u + this.wMax*w*ones(1,size(x,2)) ] - fq_;
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            q_ = this.sys0.interpolate( t );
            
            dfdx_ = [0, 0, -x(3)*cos(q_(3)) - (1 - x(3)^2/2)*sin(q_(3));...
                     0, 0, (1 - x(3)^2/2)*cos(q_(3)) + (-x(3))*sin(q_(3));...
                     0, 0, 0];
        end
        function dfdw_ = dfdw(this, x, u, w, t)
            q_ = this.sys0.interpolate( t );
             
            dfdw_ = [0; 0; this.wMax];
        end
    end
end
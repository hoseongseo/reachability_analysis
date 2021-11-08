classdef MultirotorY < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
        wMax; % maximum disturbance
        sys0; % nominal dynamics
    end
    properties (GetAccess = public, SetAccess = private)
        xN;   % nominal trajectory
        tN;   % time of xN
    end
    methods
        function this = MultirotorY(sys0, wMax)
            this = this@Dynamics.Base(sys0.Nx,sys0.Nu,sys0.Nw);
            this.wMax = wMax;
            this.sys0 = sys0;
            this.xN = sys0.xN;
            this.tN = sys0.tN;
        end
        % overloading dynamics
        function f_ = f(this, x, u, w, t)
            q_ = this.sys0.interpolate( t );
            
            fq_ = this.sys0.f(q_, u, zeros(this.Nw,1), t);
            
            cx3 = 1 - 0.5*x(3,:).^2;
%             cx3 = 1;
%             sx3 = x(3,:) - (1/6)*(x(3,:).*x(3,:).*x(3,:));
            sx3 = x(3,:);
%             cx4 = 1 - 0.5*x(4,:).^2;
%             cx4 = 1;
%             sx4 = x(4,:) - (1/6)*(x(4,:).*x(4,:).*x(4,:));
%             sx4 = x(4,:);
            
%             c3 = cx3*cos(q_(3)) - sx3*sin(q_(3));
            s3 = sx3*cos(q_(3)) + cx3*sin(q_(3));
%             c4 = cx4*cos(q_(4)) - sx4*sin(q_(4));
%             s4 = sx4*cos(q_(4)) + cx4*sin(q_(4));
            
%                 f_ = [x(4,:);...
%                     x(5,:);...
%                     x(6,:);...
%                     u(1)*c7.*s8 + this.wMax(1)*w(1) - fq_(4);...
%                     -u(1)*s7 + this.wMax(2)*w(2) - fq_(5);...
%                     u(1)*c7.*c8 - 9.8 + this.wMax(3)*w(3) - fq_(6);...
%                     Polynomial(0, zeros(size(x(1).order)));...
%                     Polynomial(0, zeros(size(x(1).order)))];
                if class(x) == "Polynomial"
                f_ = [x(2,:);...
                      -u(1)*s3 - fq_(2) + this.wMax(1)*w(1);...
                      Polynomial(u(2), zeros(size(x(1).order)));...
                      Polynomial(u(3), zeros(size(x(1).order)))];
                else
                f_ = [x(2,:);...
                      -u(1)*s3 - fq_(2) + this.wMax(1)*w(1);...
                      u(2);...
                      u(3)];
                end
%             f_ = [x(4,:);...
%                   x(5,:);...
%                   x(6,:);...
%                   u(1)*x(7,:) + this.wMax*w(1);...
%                   -u(1)*x(8,:) + this.wMax*w(2);...
%                   this.wMax*w(3)*ones(1,size(x,2));...
%                   zeros(1,size(x,2));...
%                   zeros(1,size(x,2))];
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            q_ = this.sys0.interpolate( t );
            
            cx3 = 1 - 0.5*x(3,:).^2;
%             cx3 = 1;
%             sx3 = x(3,:) - (1/6)*(x(3,:).*x(3,:).*x(3,:));
            sx3 = x(3,:);
            cx4 = 1 - 0.5*x(4,:).^2;
%             cx4 = 1;
%             sx4 = x(4,:) - (1/6)*(x(4,:).*x(4,:).*x(4,:));
            sx4 = x(4,:);
            
            c3 = cx3*cos(q_(3)) - sx3*sin(q_(3));
            s3 = sx3*cos(q_(3)) + cx3*sin(q_(3));
            c4 = cx4*cos(q_(4)) - sx4*sin(q_(4));
            s4 = sx4*cos(q_(4)) + cx4*sin(q_(4));
            
            dfdx_ = zeros(this.Nx, this.Nx);
            dfdx_(1,2) = 1.0;
            dfdx_(2,3) = -u(1)*c3;
%             dfdx_(2,4) = u(1)*c3*c4;
        end
        function dfdw_ = dfdw(this, x, u, w, t)
%             q_ = this.sys0.interpolate( t );
            dfdw_ = zeros(this.Nx, this.Nw);
            dfdw_(2,1) = this.wMax(1);
        end
    end
end
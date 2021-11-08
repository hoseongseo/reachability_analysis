classdef Multirotor < Dynamics.Base
    properties (GetAccess = private, SetAccess = private)
        wMax; % maximum disturbance
        sys0; % nominal dynamics
    end
    properties (GetAccess = public, SetAccess = private)
        xN;   % nominal trajectory
        tN;   % time of xN
    end
    methods
        function this = Multirotor(sys0, wMax)
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
            
%             cx7 = 1 - x(7,:).^2*0.5;
%             sx7 = x(7,:);
%             cx8 = 1 - x(8,:).^2*0.5;
%             sx8 = x(8,:);

            cx7 = 1 - 0.5*x(7,:).^2;
%             cx7 = 1;
%             sx7 = x(7,:) - (1/6)*(x(7,:).*x(7,:).*x(7,:));
            sx7 = x(7,:);
            cx8 = 1 - 0.5*x(8,:).^2;
%             cx8 = 1;
%             sx8 = x(8,:) - (1/6)*(x(8,:).*x(8,:).*x(8,:));
            sx8 = x(8,:);
            
            c7 = cx7*cos(q_(7)) - sx7*sin(q_(7));
            s7 = sx7*cos(q_(7)) + cx7*sin(q_(7));
            c8 = cx8*cos(q_(8)) - sx8*sin(q_(8));
            s8 = sx8*cos(q_(8)) + cx8*sin(q_(8));
            
            if class(x) == "Polynomial"
                f_ = [x(4,:);...
                    x(5,:);...
                    x(6,:);...
                    u(1)*c7.*s8  - fq_(4);...
                    -u(1)*s7 - fq_(5);...
                    u(1)*c7.*c8 - 9.8 - fq_(6);...
                    Polynomial(this.wMax(1)*w(1), zeros(size(x(1).order)));...
                    Polynomial(this.wMax(2)*w(2), zeros(size(x(1).order)))];
            else
                f_ = [x(4,:);...
                    x(5,:);...
                    x(6,:);...
                    u(1)*c7.*s8 - fq_(4) + this.wMax(1)*w(1);...
                    -u(1)*s7 - fq_(5) + this.wMax(2)*w(2);...
                    u(1)*c7.*c8 - 9.8 - fq_(6) + this.wMax(3)*w(3);...
                    zeros(1,size(x,2));...
                    zeros(1,size(x,2))];
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
            
            cx7 = 1 - 0.5*x(7,:).^2;
%             sx7 = x(7,:) - (1/6)*(x(7,:).*x(7,:).*x(7,:));
            sx7 = x(7,:);
            cx8 = 1 - 0.5*x(8,:).^2;
%             sx8 = x(8,:) - (1/6)*(x(8,:).*x(8,:).*x(8,:));
            sx8 = x(8,:);
            
            c7 = cx7*cos(q_(7)) - sx7*sin(q_(7));
            s7 = sx7*cos(q_(7)) + cx7*sin(q_(7));
            c8 = cx8*cos(q_(8)) - sx8*sin(q_(8));
            s8 = sx8*cos(q_(8)) + cx8*sin(q_(8));
            
            dfdx_ = zeros(this.Nx, this.Nx);
            dfdx_(1,4) = 1.0;
            dfdx_(2,5) = 1.0;
            dfdx_(3,6) = 1.0;
            
            dfdx_(4,7) = -u(1)*(s7.*s8);
            dfdx_(4,8) = u(1)*(c7.*c8);
            dfdx_(5,7) = -u(1)*c7;
            dfdx_(5,8) = 0.0;
            dfdx_(6,7) = -u(1)*(s7.*c8);
            dfdx_(6,8) = -u(1)*(c7.*s8);
        end
        function dfdw_ = dfdw(this, x, u, w, t)
%             q_ = this.sys0.interpolate( t );
            dfdw_ = zeros(this.Nx, this.Nw);
            dfdw_(4,1) = this.wMax(1);
            dfdw_(5,2) = this.wMax(2);
            dfdw_(6,3) = this.wMax(3);
        end
    end
end
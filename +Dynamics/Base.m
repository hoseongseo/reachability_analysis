classdef Base
    properties
        Nx;
        Nu;
        Nw;
    end
    properties (GetAccess = private, SetAccess = private)
        derivative_computed = false;
    end
    methods
        % basic constructor
        function this = Base(nx, nu, nw)
            this.Nx = nx;
            this.Nu = nu;
            this.Nw = nw;
        end
        % dynamics
        function f_ = f(this, x, u, w, t)
            f_ = zeros(this.Nx,1);
        end
        function dfdx_ = dfdx(this, x, u, w, t)
            dfdx_ = zeros(this.Nx,this.Nx);
        end
        function dfdw_ = dfdw(this, x, u, w, t)
            dfdw_ = zeros(this.Nx,this.Nw);
        end
        
        % integration
        function [x, A, D, r] = euler_forward(this, x0, u, t)
            N = length(t);
            x = zeros(this.Nx,N);
            x(:,1) = x0;
            f = zeros(this.Nx,N);
            A = zeros(this.Nx,this.Nx,N);
            D = zeros(this.Nx,this.Nw,N);
            r = zeros(this.Nx,N);
            w0 = zeros(this.Nw,1);
            
            for i = 1:N
                x_ = x(:,i);
                u_ = u(:,i);
                t_ = t(i);
                if i < N
                    dt_ = t(i+1) - t(i);
                else
                    dt_ = t(i) - t(i-1);
                end
                f_ = this.f(x_, u_, w0, t_);
                A_ = this.dfdx(x_, u_, w0, t_);
                D_ = this.dfdw(x_, u_, w0, t_);
                r_ = f_ - A_*x_;
                
                f(:,i) = f_;
                A(:,:,i) = A_;
                D(:,:,i) = D_;
                r(:,i) = r_;
                if i < N
                    x(:,i+1) = x_ + dt_*f_;
                end
                
            end
        end
        function [x, A, D, r] = rk4_forward(this, x0, u, t)
            I = eye(this.Nx);
            N = length(t);
            x = zeros(this.Nx,N);
            x(:,1) = x0;
            f = zeros(this.Nx,N);
            A = zeros(this.Nx,this.Nx,N);
            D = zeros(this.Nx,this.Nw,N);
            r = zeros(this.Nx,N);
            w0 = zeros(this.Nw,1);
            
            for i = 1:N
                u_ = u(:,i);
                t_ = t(i);
                if i < N
                    dt_ = t(i+1) - t(i);
                else
                    dt_ = t(i) - t(i-1);
                end
                
                x0_     = x(:,i);
                f0_     = this.f(x0_, u_, w0, t_);
                df0_x0  = this.dfdx(x0_, u_, w0, t_);
                df0_w   = this.dfdw(x0_, u_, w0, t_);
                
                x1_     = x0_ + 0.5*dt_*f0_;
                f1_     = this.f(x1_, u_, w0, t_);
                df1_x1  = this.dfdx(x1_, u_, w0, t_);
                df1_w   = this.dfdw(x1_, u_, w0, t_);
                df1_x0  = df1_x1 * (I + 0.5*dt_*df0_x0);
                
                x2_     = x0_ + 0.5*dt_*f1_;
                f2_     = this.f(x2_, u_, w0, t_);
                df2_x2  = this.dfdx(x2_, u_, w0, t_);
                df2_w   = this.dfdw(x2_, u_, w0, t_);
                df2_x0  = df2_x2 * (I + 0.5*dt_*df1_x0);
                
                x3_     = x0_ + dt_*f2_;
                f3_     = this.f(x3_, u_, w0, t_);
                df3_x3  = this.dfdx(x3_, u_, w0, t_);
                df3_w   = this.dfdw(x3_, u_, w0, t_);
                df3_x0  = df3_x3 * (I + dt_*df2_x0);
                
                f_ = (1/6)*(f0_ + 2*f1_ + 2*f2_ + f3_);
                A_ = (1/6)*(df0_x0 + 2*df1_x0 + 2*df2_x0 + df3_x0);
                D_ = (1/6)*(df0_w + 2*df1_w + 2*df2_w + df3_w);
                r_ = f_ - A_*x0_;
                
                f(:,i) = f_;
                A(:,:,i) = A_;
                D(:,:,i) = D_;
                r(:,i) = r_;
                if i < N
                    x(:,i+1) = x0_ + dt_*f_;
                end
            end
        end
    end
end
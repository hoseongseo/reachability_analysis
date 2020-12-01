classdef HJBequation < handle
    properties (GetAccess = public, SetAccess = private)
        grid
    end
    properties (GetAccess = private, SetAccess = private)
        sys
        wMax
        wMin
        schemeData
        HJIextraArgs
    end
    methods
        function this = HJBequation(sys, gr, hamFunc)
            %%%%
            %%%% WANNING!! dimension less than equal to 3 recommanded
            %%%%
            this.sys = sys; % dynamics
            this.grid = gr; % grid
            
            %%% prepare arguments
            this.schemeData.grid = this.grid;
            this.schemeData.uMode = 'max';
            if nargin < 3
                this.schemeData.hamFunc = @this.genericHamiltonian;
            else
                disp('Custom Hamiltonian function used!')
                this.schemeData.hamFunc = hamFunc;
            end
            this.schemeData.partialFunc = @this.genericRegulation;
            this.schemeData.accuracy = 'veryHigh';
            this.schemeData.tMode = 'forward';
            this.schemeData.system = this.sys;
            this.HJIextraArgs.visualize = false;
            this.HJIextraArgs.fig_num = 131;
            this.HJIextraArgs.deleteLastPlot = false;
            this.HJIextraArgs.quiet = false;
            
            this.wMax = ones(sys.Nw,1);
            this.wMin = -ones(sys.Nw,1);
            this.schemeData.wMax = this.wMax;
            this.schemeData.wMin = this.wMin;
        end
        
        %%%%%%%% Main function
        function V = solve(this, V0, t0, tf)
            %%%% Solve HJBequation from t0 to tf with the initial value V0
            % update disturbance bound and regulation parameter
            n = this.sys.Nx;
            statedims = repmat({':'},1,n);
            V_tmp = HJIPDE_solve(V0, [t0,tf], this.schemeData, 'none', this.HJIextraArgs);
            V = V_tmp( statedims{:}, size(V_tmp,n+1) );
        end
    end
    
    methods %(Hidden)
        %%%%%%%% Hamiltonian 
        function [hamValue, schemeData] = genericHamiltonian(this, t, data, deriv, schemeData)
            %%%% last 4 arguments are required for the solver
            nPts = prod(this.grid.shape); % number of all grid points
            n = this.sys.Nx; % number of states
            wAll = Utils.genExtrema( this.wMax ); % assuming symmetric disturbance bound
            
            % stack of state and gradiant
            x = zeros(n, nPts);
            grad = zeros(n, nPts);
            for i = 1:n
                x(i,:) = this.grid.xs{i}(:);
                grad(i,:) = deriv{i}(:);
            end
            
            %%%%% brute-force computation of hamiltonian
            % 1) compute dot product of gradient and dynamics for all possible
            % disturbances
            % 2) maximum of the dot products is the Hamiltonian
            hamVec = zeros(size(wAll,2),nPts);
            for i = 1:size(wAll,2)
                %%%%%% AT THIS TIME, INPUT IS SET AS ZERO...
                f = this.sys.f(x, zeros(this.sys.Nu,1), wAll(:,i), t);
                hamVec(i,:) = sum(f.*grad,1);
            end
            hamValue = reshape(max(hamVec,[],1), this.grid.shape);
        end
        %%%%%%%% Regulation
        function alpha = genericRegulation(this, t, data, derivMin, derivMax, schemeData, dim)
            % this function determines regularizing coefficients of viscosity term
            % alpha := max(\dot{x})
            % If this value is small, NAN would be apper. In that case, alpha must be increased enough.
            nPts = prod(this.grid.shape); % number of all grid points
            n = this.sys.Nx; % number of states
            wAll = Utils.genExtrema( this.wMax ); % assuming symmetric disturbance bound
            x = zeros(n, nPts);
            for i = 1:n
                x(i,:) = this.grid.xs{i}(:);
            end
            f = zeros(size(wAll,2)*n,nPts);
            for i = 1:size(wAll,2)
                %%%%%% AT THIS TIME, INPUT IS SET AS ZERO...
                f_ = this.sys.f(x, zeros(this.sys.Nu,1), wAll(:,i), t);
                f((i-1)*n+1:i*n,:) = f_;
            end
            tmp = max( abs(f(dim:n:end,:)),[], 1 );
            alpha = reshape(tmp, this.grid.shape );
        end
    end
end
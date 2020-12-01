classdef FunnelLTV < handle
    properties
        Q
        t
        Q_basis
        a_basis
    end
    properties (GetAccess = private, SetAccess = private)
        sys
        composition_method = "determinant";
        internal_composition = false;
        eps = 1e-10;
    end
    methods
        function this = FunnelLTV(sys, Q0, tk, arg1, arg2)
            
            if nargin == 4
                this.composition_method = arg1;
            elseif nargin == 5
                this.composition_method = arg1;
                this.internal_composition = arg2;
            end
            
            nx = sys.Nx;
            nw = sys.Nw;
            N = length(tk);
            
            this.t = tk;
            this.Q = zeros(nx,nx,N);
            if this.composition_method == "determinant"
                this.Q(:,:,1) = Composition.MVOE(Q0);
            elseif this.composition_method == "trace"
                this.Q(:,:,1) = Composition.MTOE(Q0);
            end
            this.Q_basis = cell(1,N);
            this.Q_basis{1} = Q0;
            this.a_basis = cell(1,N);
            this.a_basis{1} = 1.0;
             
            this.sys = sys;
            
            % fourth dimension is now disturbance channel
            Nk = zeros(nx,nx,N-1,nw);
            
            P = eye(nx); % this should be initialized outside for consecutive propagation
            for k = 1:N-1
                t_ = tk(k);
                dt_ = tk(k+1) - tk(k);

                A_ = sys.dfdx(zeros(sys.Nx,1), zeros(sys.Nu,1), zeros(sys.Nw,1), t_);
                D_ = sys.dfdw(zeros(sys.Nx,1), zeros(sys.Nu,1), zeros(sys.Nw,1), t_);
                TF_ = expm(A_*dt_);
                
                %%% sets due to disturbances
                rank_deficient = rank(A_) < this.sys.Nx;
                if rank_deficient
                    disp(['A is singular @ t=',num2str(tk(k)),' (k=',num2str(k),')'])
                end
                for j = 1:nw
                    M_ = (D_(:,j)*D_(:,j)');
                    if rank_deficient
                        N_ = (M_ + expm(-A_*dt_)*M_*expm(-A_'*dt_))*0.5*dt_; % approximation
                    else
                        N_ = lyap(A_, expm(-A_*dt_)*M_*expm(-A_'*dt_) - M_); % true lyapunov
                    end
                    Nk(:,:,k,j) = dt_*(N_ + this.eps*dt_*(P*P'));
                end
                P = expm(A_*dt_)*P;
                
                %%% prepare for basis
                N_basis_prev = size(this.Q_basis{k},3); % start from the previous basis
                N_basis = N_basis_prev; 
                if this.internal_composition
                    % the number of new ellipsoid is always "1"
                    N_basis = N_basis + 1;
                else
                    % the number of new ellipsoid is always "nw"
                    N_basis = N_basis + nw;
                end
                Qtmp = zeros(nx,nx,N_basis);
                Qtmp(:,:,1:N_basis_prev) = this.Q_basis{k}; % get the previous one
                
                % effect of disturbances
                Ntmp = squeeze(Nk(:,:,k,:));
                if this.internal_composition
                    % compose disturbance ellipsoids at k-th index(slight conservative)
                    if this.composition_method == "determinant"
                        Qsub = Composition.MVOE(Ntmp);
                    elseif this.composition_method == "trace"
                        Qsub = Composition.MTOE(Ntmp);
                    end
                    Qtmp(:,:,N_basis_prev+1) = Qsub;
                else
                    % just stack disturbance ellipsoids (what I intended)
                    Qtmp(:,:,N_basis_prev+(1:nw)) = Ntmp;
                end
                
                % this would be right position to transform
                for i = 1:size(Qtmp,3)
                    Qtmp(:,:,i) = TF_ * Qtmp(:,:,i) * TF_';
                end
                
                % compose ellipsoids now and use for next step
                this.Q_basis{k+1} = Qtmp;
                if this.composition_method == "determinant"
                    [this.Q(:,:,k+1), a_k] = Composition.MVOE( this.Q_basis{k+1} );
                elseif this.composition_method == "trace"
                    [this.Q(:,:,k+1), a_k] = Composition.MTOE( this.Q_basis{k+1} );
                elseif this.composition_method == "sum"
                    this.Q(:,:,k+1) = sum( this.Q_basis{k+1}, 3 );
                    a_k = ones(1,size(this.Q_basis{k+1}, 3));
                end
                this.a_basis{k+1} = a_k;
            end
            
        end
    end
end
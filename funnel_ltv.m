function Q = funnel_ltv( sys, tk, Q0, args )
% sys: dynamics
% tk: time horizon of interest (discretized)
% Q0: shape matrix for the initial set
% args: extra arguments detailed as
if nargin < 4
    args.internal_composition = false;
    args.composition_method = "MVOE";
    args.eps = 1e-10;
end

N = length(tk);

Q = zeros(sys.Nx,sys.Nx,N);
Q(:,:,1) = Q0;
Q_basis = cell(1,N);
Q_basis{1} = Q0;
a_basis = cell(1,N);
a_basis{1} = 1.0;

% the fourth dimension represents disturbance channel
Nk = zeros(sys.Nx,sys.Nx,N-1,sys.Nw);

% state-transition matrix
P = eye(sys.Nx); % this should be initialized outside for consecutive propagation
for k = 1:N-1
    t_ = tk(k);
    dt_ = tk(k+1) - tk(k);
    
    A_ = sys.dfdx(zeros(sys.Nx,1), zeros(sys.Nu,1), zeros(sys.Nw,1), t_);
    D_ = sys.dfdw(zeros(sys.Nx,1), zeros(sys.Nu,1), zeros(sys.Nw,1), t_);
    TF_ = expm(A_*dt_);
    
    %%% sets due to disturbances
    rank_deficient = rank(A_) < sys.Nx;
    if rank_deficient
        disp(['A is singular @ t=',num2str(tk(k)),' (k=',num2str(k),')'])
    end
    for j = 1:sys.Nw
        M_ = (D_(:,j)*D_(:,j)');
        if rank_deficient
            N_ = (M_ + expm(-A_*dt_)*M_*expm(-A_'*dt_))*0.5*dt_; % approximation
        else
            N_ = lyap(A_, expm(-A_*dt_)*M_*expm(-A_'*dt_) - M_); % true lyapunov
        end
        Nk(:,:,k,j) = dt_*(N_ + args.eps*dt_*(P*P'));
    end
    P = expm(A_*dt_)*P;
    
    %%% prepare for basis
    N_basis_prev = size(Q_basis{k},3); % start from the previous basis
    N_basis = N_basis_prev;
    if args.internal_composition
        % the number of new ellipsoid is always "1"
        N_basis = N_basis + 1;
    else
        % the number of new ellipsoid is always "nw"
        N_basis = N_basis + sys.Nw;
    end
    Qtmp = zeros(sys.Nx,sys.Nx,N_basis);
    Qtmp(:,:,1:N_basis_prev) = Q_basis{k}; % get the previous one
    
    % effect of disturbances
    Ntmp = squeeze(Nk(:,:,k,:));
    if args.internal_composition
        % compose disturbance ellipsoids at k-th index(slight conservative)
        if args.composition_method == "MVOE"
            Qsub = Composition.MVOE(Ntmp);
        elseif args.composition_method == "MTOE"
            Qsub = Composition.MTOE(Ntmp);
        end
        Qtmp(:,:,N_basis_prev+1) = Qsub;
    else
        % just stack disturbance ellipsoids
        Qtmp(:,:,N_basis_prev+(1:sys.Nw)) = Ntmp;
    end
    
    % this would be right position to transform
    for i = 1:size(Qtmp,3)
        Qtmp(:,:,i) = TF_ * Qtmp(:,:,i) * TF_';
    end
    
    % compose ellipsoids now and use for next step
    Q_basis{k+1} = Qtmp;
    if args.composition_method == "MVOE"
        [Q(:,:,k+1), a_k] = Composition.MVOE( Q_basis{k+1} );
    elseif args.composition_method == "MTOE"
        [Q(:,:,k+1), a_k] = Composition.MTOE( Q_basis{k+1} );
    end
    a_basis{k+1} = a_k;
end
end
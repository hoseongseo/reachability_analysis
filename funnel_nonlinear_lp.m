function [region, cost_hist, rate_hist] = funnel_nonlinear_lp(sys, tk, Q0, N, args)
% sys: dynamics
% tk: time horizon of interest (discretized)
% Q0: shape matrix for the initial set
% N: target order of the approximated value function
% args: extra arguments detailed as
if nargin < 5
    args.max_iter = 10;
    args.ftol = 1e-5;
    args.plot_cost = false;
end

% declare state and disturbance
xp(sys.Nx,1) = Polynomial;
for i = 1:sys.Nx
    ord_ = zeros(1,sys.Nx+sys.Nw);
    ord_(i) = 1;
    xp(i) = Polynomial(1, ord_);
end
wp(sys.Nw,1) = Polynomial;
for i = 1:sys.Nw
    ord_ = zeros(1,sys.Nx+sys.Nw);
    ord_(sys.Nx+i) = 1;
    wp(i) = Polynomial(1, ord_);
end

% domain of interest (initialized by the funnel of LTV system)
Q_ltv = funnel_ltv(sys, tk, Q0);
region = struct('q', zeros(sys.Nx,length(tk)), 'Q', Q_ltv);

% degrees of the approximated value functions
poly_degree = Polynomial.integer_grid(N);
quad_degree = poly_degree(sum(poly_degree, 2) <= 2, :);
quad_idx_1st_order = sum(quad_degree,2) == 1;
quad_idx_2nd_order = sum(quad_degree,2) == 2;
quad_idx_skew = max(quad_degree,[],2) == 1;

% initial coefficients
V0 = xp'*inv(Q0)*xp - 1;
c0 = Polynomial.expand_matrix(N, V0.order(:,1:sys.Nx))*V0.coeff;
c_tilde0 = c0(sum(poly_degree, 2) <= 2);

% order of dynamics
fp = sys.f(xp, zeros(sys.Nu,1), wp, tk(1));
N_dyn = max(reshape([fp.order], [numel([fp.order])/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);
Nx = N_dyn(1:sys.Nx); % order in state
Nw = N_dyn(sys.Nx + (1:sys.Nw)); % order in disturbance

% precomputation of matrices
D_N = Polynomial.derivative_matrix( N );
invB_Nw = inv(Polynomial.bernstein_transform_matrix( Nw ));
B_N = Polynomial.bernstein_transform_matrix( N );

P_N_Nx = Polynomial.product_matrix(N, Nx);
I_N = eye(prod(N+1));
C_N_Nx = Polynomial.expand_matrix( N+Nx, Polynomial.integer_grid(N));
E_N = Polynomial.expand_matrix( N, quad_degree );

B_NNx = (Polynomial.bernstein_transform_matrix( N+Nx ));
h_N = prod((1./(poly_degree+1)).*( ones(size(poly_degree)).^(poly_degree+1) - (-ones(size(poly_degree))).^(poly_degree+1) ),2);

% Hamiltonian-bounding matrices over time
H = zeros(prod(N+Nx+1), prod(N+1), prod(Nw+1), length(tk)-1);
for k = 1:length(tk)-1
    f_ = sys.f(xp, zeros(sys.Nu,1), wp, tk(k));
    S_ = Polynomial.dyn_coeff( f_, sys.Nx, sys.Nw, N_dyn );
    
    mat = zeros(prod(N+Nx+1), prod(N+1), prod(Nw+1));
    for m = 1:prod(Nw+1)
        for j = 1:sys.Nx
            mat(:,:,m) = mat(:,:,m) + P_N_Nx * kron(I_N, S_(:,m,j)) * D_N(:,:,j);
        end
    end
    
    for m = 1:prod(Nw+1)
        for i = 1:prod(Nw+1)
            H(:,:,m,k) = H(:,:,m,k) + invB_Nw(m,i) * mat(:,:,i);
        end
    end
end

% volume of the initial guess
vol = 0.0;
for k = 2:length(tk)
    vol = vol + log(det( region.Q(:,:,k) ));
end
cost_hist = vol;
rate_hist = [];

% MAIN LOOP
for iter = 1:args.max_iter
    %% Domain transform matrices over time
    % current domain of interest
    region_ = region(end);
    
    T_NNx = zeros( prod(N+Nx+1), prod(N+Nx+1), length(tk) );
    T_N = zeros( prod(N+1), prod(N+1), length(tk) );
    for k = 1:length(tk)
        sqrtQ_ = region_.Q(:,:,k)^(1/2);
        a_ = sqrt(sum(sqrtQ_.^2,2));
        b_ = region_.q(:,k);
        T_NNx(:,:,k) = Polynomial.domain_transform_matrix( N+Nx, a_, b_ );
        T_N(:,:,k) = Polynomial.domain_transform_matrix( N, a_, b_ );
    end
    
    %% Solve LP
    cvx_begin
        variable c( prod(N+1), length(tk) )
        variable c_tilde( (sys.Nx+1)*(sys.Nx+2)/2, length(tk) )
    
        cost = 0.0;
        for k = 2:length(tk)
            cost = cost - h_N'*(T_N(:,:,k)*E_N*c_tilde(:,k));
        end
        minimize cost
    
        subject to
            c(:,1) == c0
            c_tilde(:,1) == c_tilde0
            for k = 1:length(tk)-1
                for m = 1:prod(Nw+1)
                    B_NNx\T_NNx(:,:,k)*( C_N_Nx*(c(:,k+1)-c(:,k)) + (tk(k+1)-tk(k))*H(:,:,m,k)*c(:,k) ) <= 0
                end
                B_N\T_N(:,:,k+1)*(E_N*c_tilde(:,k+1) - c(:,k+1)) <= 0
            end
    cvx_end
    
    %% Post processes
    % extract ellipsoid from c_tilde
    c_tilde_conv = c_tilde.*repmat( ones((sys.Nx+1)*(sys.Nx+2)/2,1) - 0.5*quad_idx_skew, [1,length(tk)] );
    r = c_tilde_conv(1,:);
    d = c_tilde_conv(quad_idx_1st_order, :);
    vechA = c_tilde_conv(quad_idx_2nd_order, :);
    
    % domain to be updated
    region_new = region_;
    for k = 1:length(tk)
        r_ = r(k);
        d_ = d(:,k);
        A_ = reshape(Utils.duplication_matrix(sys.Nx)*vechA(:,k), [sys.Nx,sys.Nx]);
        
        q_ = -A_\d_;
        Q_ = (d_'*(A_\d_)-r_)*inv(A_);
        
        region_new.q(:,k) = q_;
        region_new.Q(:,:,k) = Q_;
    end
    
    % current cost
    vol = 0.0; 
    for k = 2:length(tk)
        vol = vol + log(det( region_new.Q(:,:,k) ));
    end
    dec_rate = (cost_hist(end) - vol) / cost_hist(end);
    disp(['Iteration #', num2str(iter),...
        ': cost = ', num2str(vol),...
        ', dcost = ', num2str(cost_hist(end)-vol),...
        ', rate = ', num2str(dec_rate)]);
    
    cost_hist = [cost_hist, vol];
    rate_hist = [rate_hist, dec_rate];
    
    % visualization
    if args.plot_cost
        figure(123)
        subplot(2,1,1)
        title('$\textbf{Convergence characteristics of SLP}$')
        cla; hold on; grid on;
        plot(cost_hist,'b*-')
        axis tight;
        ylabel('$cost$')
        
        subplot(2,1,2)
        cla; hold on; grid on;
        plot(rate_hist,'b*-')
        axis tight;
        ylabel('$rate$')
        xlabel('Iterations')
        drawnow
    end
    
    %% Convergence check & domain update
    if abs(dec_rate) < args.ftol
        disp(['Converged (rate = ',num2str(dec_rate),' / ftol = ',num2str(args.ftol),')'])
        break;
    end
    region(iter+1) = region_new;
end
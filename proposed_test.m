clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))

% given initial set
Q0 = diag([0.05; 0.05])^2;

% disturbance bound
wMax = 0.1;

% sphere
S = Utils.Sphere(1,200);

%% system
t = linspace(0,1,101);
% nominal dynamics with the initial condition
sys0 = Dynamics.LotkaVolterraNominal([1.2; 1.1], t);
% nonlinear system shifted to the origin
sys = Dynamics.LotkaVolterra(sys0, wMax);


%% HJB equation
load('V_lotka_volterra.mat')
X = cell(1,length(t));
gr = createGrid([-0.06, -0.06], [0.06, 0.06], [201,201]);
for k = 1:length(t)
    X{k} = getLevelSet( gr, V(:,:,k), 0.0 );
end

%% Polynomial variable
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


%% Initialization
% utilize LTV version
funnel_ltv = FunnelLTV(sys, Q0, t);

% domain of interest
R = struct('q', zeros(sys.Nx,length(t)), 'Q', zeros(sys.Nx,sys.Nx,length(t)));
R.Q = funnel_ltv.Q;

% target order of value function
N = [3,3]; 
poly_degree = integer_grid(N);

% % quadradic function
quad_degree = poly_degree(sum(poly_degree, 2) <= 2, :);

% initial coefficient
V0 = xp'*inv(Q0)*xp - 1;
c0 = Polynomial.expand_matrix(N, V0.order(:,1:sys.Nx))*V0.coeff;
c_tilde0 = c0(sum(poly_degree, 2) <= 2);

% order of dynamics
fp = sys.f(xp, zeros(sys.Nu,1), wp, t(1));
N_dyn = max(reshape([fp.order], [numel([fp.order])/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);
Nx = N_dyn(1:sys.Nx); % maximum order of state
Nw = N_dyn(sys.Nx + (1:sys.Nw)); % maximum order of disturbance

% preparation
D_N = Polynomial.derivative_matrix( N );
B_Nw = inv(Polynomial.bernstein_transform_matrix( Nw ));
B_N = Polynomial.bernstein_transform_matrix( N );

P_N_Nx = Polynomial.product_matrix(N, Nx);
I_N = eye(prod(N+1));
C_N_Nx = Polynomial.expand_matrix( N+Nx, integer_grid(N));
E_N = Polynomial.expand_matrix( N, quad_degree );

B_NNx = (Polynomial.bernstein_transform_matrix( N+Nx ));
h_N = prod((1./(poly_degree+1)).*( ones(size(poly_degree)).^(poly_degree+1) - (-ones(size(poly_degree))).^(poly_degree+1) ),2);

quad_idx_1st_order = sum(quad_degree,2) == 1;
quad_idx_2nd_order = sum(quad_degree,2) == 2;
quad_idx_skew = max(quad_degree,[],2) == 1;

% Hamiltonian bounding matrix over time 
H = zeros(prod(N+Nx+1), prod(N+1), prod(Nw+1), length(t)-1);
for k = 1:length(t)-1
    f_ = sys.f(xp, zeros(sys.Nu,1), wp, t(k));
    S_ = Polynomial.dyn_coeff( f_, sys.Nx, sys.Nw, N_dyn );
    
    mat = zeros(prod(N+Nx+1), prod(N+1), prod(Nw+1));
    for m = 1:prod(Nw+1)
        for j = 1:sys.Nx
            mat(:,:,m) = mat(:,:,m) + P_N_Nx * kron(I_N, S_(:,m,j)) * D_N(:,:,j);
        end
    end
    
    for m = 1:prod(Nw+1)
        for i = 1:prod(Nw+1)
            H(:,:,m,k) = H(:,:,m,k) + B_Nw(m,i) * mat(:,:,i);
        end
    end
end

% volume of the initial guess
vol = 0.0;
for k = 2:length(t)
    vol = vol + log(det( R.Q(:,:,k) ));
end

costHist = vol;
rateHist = [];

ftol = 1e-5;

for iter = 1:10
%% Domain transform matrix
T_NNx = zeros( prod(N+Nx+1), prod(N+Nx+1), length(t) );
T_N = zeros( prod(N+1), prod(N+1), length(t) );
for k = 1:length(t)
    sqrtQ_ = R.Q(:,:,k)^(1/2);
    a_ = sqrt(sum(sqrtQ_.^2,2));
    b_ = R.q(:,k);
    T_NNx(:,:,k) = Polynomial.domain_transform_matrix( N+Nx, a_, b_ );
    T_N(:,:,k) = Polynomial.domain_transform_matrix( N, a_, b_ );
end

%% Solve for all
cvx_begin
    variable c( prod(N+1), length(t) )
    variable c_tilde( (sys.Nx+1)*(sys.Nx+2)/2, length(t) )

    cost = 0.0;
    for k = 2:length(t)
        cost = cost - h_N'*(T_N(:,:,k)*E_N*c_tilde(:,k));
    end
    minimize cost
    
    subject to
        c(:,1) == c0
        c_tilde(:,1) == c_tilde0
        for k = 1:length(t)-1
            for m = 1:prod(Nw+1)
                B_NNx\T_NNx(:,:,k)*( C_N_Nx*(c(:,k+1)-c(:,k)) + (t(k+1)-t(k))*H(:,:,m,k)*c(:,k) ) <= 0
            end
            B_N\T_N(:,:,k+1)*(E_N*c_tilde(:,k+1) - c(:,k+1)) <= 0
        end
    
cvx_end

% %% Polynomial approximation (Step 1)
% cvx_begin
%     variable c( prod(N+1), length(t) )
%     
%     cost1 = 0.0;
%     for k = 2:length(t)
%         cost1 = cost1 - h_N'*(T_N(:,:,k)*c(:,k));
%     end
%     minimize cost1
% 
%     subject to
%         c(:,1) == c0
%         for k = 1:length(t)-1
%             for m = 1:prod(Nw+1)
%                 B_NNx\T_NNx(:,:,k)*( C_N_Nx*(c(:,k+1)-c(:,k)) + (t(k+1)-t(k))*H(:,:,m,k)*c(:,k) ) <= 0
%             end
%         end
% cvx_end
% 
% %% Quadratic approximation (Step 2)
% cvx_begin
%     variable c_tilde( (sys.Nx+1)*(sys.Nx+2)/2, length(t) )
%     
%     cost2 = 0.0;
%     for k = 2:length(t)
%         cost2 = cost2 - h_N'*(T_N(:,:,k)*E_N*c_tilde(:,k));
%     end
%     minimize cost2
%     
%     subject to
%         c_tilde(:,1) == c_tilde0
%         for k = 2:length(t)
%             B_N\T_N(:,:,k)*(E_N*c_tilde(:,k) - c(:,k)) <= 0
%         end
% cvx_end

%% post process: extract ellipsoid from c_tilde
c_tilde_conv = c_tilde.*repmat( ones((sys.Nx+1)*(sys.Nx+2)/2,1) - 0.5*quad_idx_skew, [1,length(t)] );
r = c_tilde_conv(1,:);
d = c_tilde_conv(quad_idx_1st_order, :);
vechA = c_tilde_conv(quad_idx_2nd_order, :);

R_new = R;
for k = 1:length(t)
    r_ = r(k);
    d_ = d(:,k);
    A_ = reshape(duplication_matrix(sys.Nx)*vechA(:,k), [sys.Nx,sys.Nx]);
    
    q_ = -A_\d_;
    Q_ = (d_'*(A_\d_)-r_)*inv(A_);
    
    R_new.q(:,k) = q_;
    R_new.Q(:,:,k) = Q_;
end

vol = 0.0; % current cost
for k = 2:length(t)
    vol = vol + log(det( R_new.Q(:,:,k) ));
end
dec_rate = (costHist(end) - vol) / costHist(end);
disp(['Iteration #', num2str(iter),...
    ': cost = ', num2str(vol),...
    ', dcost = ', num2str(costHist(end)-vol),...
    ', rate = ', num2str(dec_rate)]);

costHist = [costHist, vol];
rateHist = [rateHist, dec_rate];


%% visualization
figure(331)
cla; hold on; grid on; axis tight;
for k = round(linspace(1,length(t),11))
    tmp = sys.xN(:,k) + R_new.q(:,k) + R_new.Q(:,:,k)^(1/2) * S.x;
    plot(tmp(1,:), tmp(2,:), 'r', 'linewidth', 2)
    
%     tmp = sys.xN(:,k) + X{k};
%     plot(tmp(1,:), tmp(2,:), 'k--', 'linewidth', 2)
end
drawnow

figure(34)
subplot(2,1,1)
cla; hold on; grid on;
plot(costHist,'b*-')
axis tight;
ylabel('$cost$')
subplot(2,1,2)
cla; hold on; grid on;
plot(rateHist,'b*-')
axis tight;
ylabel('$rate$')
xlabel('Iteration')
drawnow

%% Convergence check & Domain update
if abs(dec_rate) < ftol
    res = 1; % converged
    disp(['Converged (rate = ',num2str(dec_rate),' / ftol = ',num2str(ftol),')'])
    break;
end

R = R_new;
end

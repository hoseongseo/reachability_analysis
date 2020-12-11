clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

% given initial set
Q0 = diag([0.25, 0.25]);

% Q0 = diag([0.05; 0.05])^2;

% LTV system
t = linspace(0.0,1.0,31);
A = @(t) [-0.8*t + 0.5, cos(1.5*t + 2); 0.5*t^(2/3), -2*exp(-0.7*t)];
D = @(t) [0.4*cos(t), -0.4*t^2; 0.08*t,  2.8*cos(3*t)];
sys = Dynamics.LTV(2, 0, 2, A, D); 

Q_ltv = funnel_ltv(sys, t, Q0);

% wMax = 0.1;
% sys0 = Dynamics.LotkaVolterraNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
% sys = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin


%%%% HJB equation
% grid for HJB equation
nGrid = [201, 201];
minGrid = [-3, -3];
maxGrid = [3, 3];
gr = createGrid(minGrid, maxGrid, nGrid);
V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) - 1;
X0 = Utils.get_level_set(gr, V0, 0.0);

% nGrid = [201, 201];
% minGrid = [-0.1, -0.1];
% maxGrid = [0.1, 0.1];
% gr = createGrid(minGrid, maxGrid, nGrid);
% V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) - 1;
% X0 = Utils.get_level_set(gr, V0, 0.0);
%%
% solve
hjb_equation = HJBequation(sys, gr);
V = zeros([size(V0), length(t)]);
V(:,:,1) = V0;
for i = 1:length(t)-1
    V(:,:,i+1) = hjb_equation.solve(V(:,:,i), t(i), t(i+1));
end

% extract zero-level set
X = cell(1,length(t));
X{1} = X0;
for i = 2:length(t)
    X{i} = Utils.get_level_set(gr, V(:,:,i), 0.0);
end

%% what can be the input?
% 1) order of value function: N
% 2) order of slack polynomial: N_lamba


%% polynomial variables
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

fp = sys.f(xp, zeros(sys.Nu,1), wp, t(1));
N_f = max(reshape([fp.order], [numel([fp.order])/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);
N_lambda = [2,2,1,1];
N_s = [2,2,1,1];
N = [2,2,zeros(1,sys.Nw)];

% order of dV*f
% N_dVf = N + N_f;

%%
% order of value*dynamics or slack polynomial determines the maximum degree
N_max = max(max(N+N_f, N+N_lambda), N_s);

E_N = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N));
E_N_N_lambda = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N + N_lambda));
E_N_N_f = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N + N_f));
E_N_s = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_s));
% E_N_lambda2 = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a+N_lambda));

B_N_max = Polynomial.bernstein_transform_matrix(N_max);
B_N_lambda = Polynomial.bernstein_transform_matrix(N_lambda);
B_N_s = Polynomial.bernstein_transform_matrix(N_s);

% Hamiltonian-bounding matrices over time
D = Polynomial.derivative_matrix(N);
P_N_N_f = Polynomial.product_matrix(N, N_f);
P_N_N_lambda = Polynomial.product_matrix(N, N_lambda);

M = zeros( prod(N+N_f+1), prod(N+1), length(t)-1 );
for k = 1:length(t)-1
    mat = zeros(prod(N+N_f+1), prod(N+1));
    f_ = sys.f(xp, zeros(sys.Nu,1), wp, t(k));
    for j = 1:sys.Nx
        mat = mat + P_N_N_f * kron( eye(prod(N+1)), Polynomial.expand_matrix(N_f, f_(j).order)*f_(j).coeff ) * D(:,:,j);
    end
    M(:,:,k) = mat;
end

%%
% initial coeffs
V0 = xp'*inv(Q0)*xp - 1;
c0 = Polynomial.expand_matrix(N, V0.order(:,1:sys.Nx))*V0.coeff;

% %% coefficient dynamics
% c = sym('c', [prod(N+1),1], 'real');
% lambda = sym('lambda', [prod(N_lambda+1),1], 'real');
%% domain transform
T_N_max = zeros(prod(N_max+1), prod(N_max+1), length(t));
T_N_lambda = zeros(prod(N_lambda+1), prod(N_lambda+1), length(t));
T_N_s = zeros(prod(N_s+1), prod(N_s+1), length(t));
for k = 1:length(t)
    sqrtQ_ = Q_ltv(:,:,k)^(1/2);
    a_ = sqrt(sum(sqrtQ_.^2,2));
    b_ = zeros(sys.Nx,1);
    T_N_max(:,:,k) = Polynomial.domain_transform_matrix( N_max, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
    T_N_lambda(:,:,k) = Polynomial.domain_transform_matrix( N_lambda, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
    T_N_s(:,:,k) = Polynomial.domain_transform_matrix( N_s, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
end

%%
CUDPArgs.dt = t(2)-t(1);
CUDPArgs.M = M;
CUDPArgs.E_N = E_N;
CUDPArgs.E_N_N_f = E_N_N_f;
CUDPArgs.E_N_N_lambda = E_N_N_lambda;
CUDPArgs.E_N_s = E_N_s;
CUDPArgs.rho_vol = 1;
CUDPArgs.rho_lambda = 1;
CUDPArgs.rho_s = 1;
CUDPArgs.B_N_s = B_N_s;
CUDPArgs.B_N_lambda = B_N_lambda;
CUDPArgs.B_N_max = B_N_max;
CUDPArgs.T_N_s = T_N_s;
CUDPArgs.T_N_lambda = T_N_lambda;
CUDPArgs.T_N_max = T_N_max;
CUDPArgs.P_N_N_lambda = P_N_N_lambda;
CUDPArgs.N = N;
CUDPArgs.N_s = N_s;
CUDPArgs.N_lambda = N_lambda;
CUDPArgs.c0 = c0;

CUDPOptions = struct(...
    'rho', 1e-5,...
    'drho', 1.0,...
    'rhoFactor', 1.6,...
    'rhoMax', 1e+10,...
    'rhoMin', 1e-6,...
    'tolGrads', 10.^linspace(-2,-5,7),...
    'tolCosts', 10.^linspace(-2,-5,7),...
    'tolConsts', 10.^linspace(-2,-5,7),...
    'alphas', 10.^linspace(0,-3,5),...
    'maxIter', 1000,...
    'mu', 0.2,...
    'lambda', 0.005,...
    'phi', 0.01,...
    'verbosity', 1);
dyn = @(x,u,idx) Problem.ValueCoeffs.dynamics(x,u,idx,CUDPArgs);
cost = @(x,u,idx) Problem.ValueCoeffs.cost(x,u,idx,CUDPArgs);
const = @(x,u,idx) Problem.ValueCoeffs.constraint(x,u,idx,CUDPArgs);
dyncst = struct('dynamics', dyn, 'cost', cost, 'constraint', const);

%%

c_tilde0 = (B_N_max\T_N_max(:,:,1))*E_N * c0;
u0 = [repmat(c0, [1,length(t)-1]); zeros(prod(N_lambda+1), length(t)-1); 0.01*ones(prod(N_s+1), length(t)-1)];

% prob = Planner.CUDP(dyncst, c0, zeros(prod(N_lambda+1),length(t)-1), CUDPOptions);
prob = Planner.CUDP(dyncst, c_tilde0, u0, CUDPOptions);
prob.solve

%%
c_from_x = zeros(prod(N+1),length(t));
c_from_u = zeros(prod(N+1),length(t));
X_cudp1 = cell(1,length(t));
X_cudp2 = cell(1,length(t));
for k = 1:length(t)
%     V_ = Polynomial(prob.xN(:,k), Polynomial.integer_grid(N));
    c_tilde_ = prob.xN(:,k);
    c_ = pinv( (B_N_max\T_N_max(:,:,k))*E_N )*c_tilde_;
    c_from_x(:,k) = c_;
    if k < length(t)
        c_from_u_ = prob.uN(1:prod(N+1),k);
        c_from_u(:,k) = c_from_u_;
        lambda_ = prob.uN(prod(N+1) + (1:prod(N_lambda+1)), k);
        s_ = prob.uN(prod(N+1) + prod(N_lambda+1) + (1:prod(N_s+1)), k);
        
        V_ = Polynomial(c_from_u_, Polynomial.integer_grid(N(:,1:sys.Nx)));
        X_ = Utils.get_level_set( gr, V_.eval(gr), 0.0 );
        X_cudp1{k} = X_;
    else
        X_cudp1{k} = X_cudp1{k-1};
    end
    V_ = Polynomial(c_, Polynomial.integer_grid(N(:,1:sys.Nx)));
    X_ = Utils.get_level_set( gr, V_.eval(gr), 0.0 );
    X_cudp2{k} = X_;
end


figure(11)
cla; hold on; grid on;
for k = round(linspace(1,length(t),5))
    %             tmp = Q_ltv(:,:,k)^(1/2) * Utils.Sphere(1,200).x;
    %             plot3(tmp(1,:), tk(k)*ones(1,size(tmp,2)), tmp(2,:), 'g', 'linewidth', 2)
    
    tmp = X{k};
    plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'k', 'linewidth', 2);
    
%     tmp = X_poly{k};
%     plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'r', 'linewidth', 2)
    
    tmp = X_cudp1{k}; % c_from_u
    plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'b', 'linewidth', 2)
    
    tmp = X_cudp2{k}; % c_from_x
    plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'r--', 'linewidth', 2)
end
xlabel('$x_1$')
zlabel('$x_2$')
ylabel('$t$ [s]')
view([128,24])
% view([-180,0])

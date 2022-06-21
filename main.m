clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS')) % HJB equation solver
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

set(groot, 'DefaultFigureColor', 'w')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
set(groot, 'DefaultAxesFontSize', 16)

% original dynamics
f = @(x,w)...
    [(w+3).*x(1,:).*(1-x(2,:));...
    (w+3).*x(2,:).*(x(1,:)-1)];
dfdx = @(x,w)...
    [(w+3)*(1-x(2)), -(w+3)*x(1);...
    (w+3)*x(2), (w+3)*(x(1)-1)];
dfdw = @(x,w)...
    [x(1)*(1-x(2));...
    x(2)*(x(1)-1)];

q = [1.2; 1.1];
Q = 0.05^2*eye(2);
wMax = 0.2;

% polynomial expression of dynamics
xp = [Polynomial(1,[1,0,0]); Polynomial(1,[0,1,0])];
wp = Polynomial(1,[0,0,1]);
fp = f(xp, wp);
t = linspace(0,1,101);

% dynamics
sys0 = Dynamics.LotkaVolterraNominal(q, t); % nominal dynamics with the initial condition
sys = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin

%% Level-set method
nGrid = [401, 401];
minGrid = [-0.1, -0.1];
maxGrid = [0.1, 0.1];
gr = createGrid(minGrid, maxGrid, nGrid);
V0 = gr.xs{1}.*gr.xs{1}/Q(1,1) + gr.xs{2}.*gr.xs{2}/Q(2,2) - 1;
X0 = Utils.get_level_set(gr, V0, 0.0);

% solve HJB PDE
hjb_equation = HJBequation(sys, gr);
V = zeros([size(V0), length(t)]);
V(:,:,1) = V0;

tic
for i = 1:length(t)-1
    V(:,:,i+1) = hjb_equation.solve(V(:,:,i), t(i), t(i+1));
end
ctime_levelset = toc;

% extract zero-level set
X = cell(1,length(t));
X{1} = X0;
for i = 2:length(t)
    X{i} = Utils.get_level_set(gr, V(:,:,i), 0.0);
end

%% Linearization-based method
W = wMax^2;

Q_lin = zeros(2,2,length(t));
Q_lin(:,:,1) = Q;
H_lin = zeros(2,length(t));
I = eye(2);
tic
for k = 1:length(t)-1
    q_ = sys.xN(:,k);
    Q_ = Q_lin(:,:,k);
    H_ = H_lin(:,k);
    
    dt_ = t(k+1)-t(k);
    A0_ = dfdx(q_, 0);
    D0_ = dfdw(q_, 0);
    f0_ = f(q_, 0);
    
    Ad_ = expm(A0_*dt_);
    Dd_ = A0_\(Ad_-I)*D0_;

    Q_lin(:,:,k+1) = Ad_*Q_*Ad_' + Dd_*H_'*Ad_' + Ad_*H_*Dd_' + Dd_*W*Dd_';
    H_lin(:,k+1) = Ad_*H_ + Dd_*W; 
end
ctime_linearization = toc;

%% SOS-program-based method
args.max_iter = 10;
args.rho = 3.0;
args.ftol = 5e-4;
args.plot_cost = false;
tic;
[res_sos, cost_sos, rate_sos] = sos_program(sys, t, Q, args);
ctime_sos = toc;

Q_sos = res_sos(end).step2;

%% Nonlinear-optimization-based method
tic
Q_nonlin = nonlinear_optimization(sys, t, Q, wMax);
ctime_nonlinopt = toc;

%% Stochastic propagation
tic
P_gmm = stochastic_propagation(q, 0.0005*eye(2), wMax, t);
ctime_stochastic = toc;

%% Proposed
Q_proposed = zeros(2,2,length(t));
Q_proposed(:,:,1) = Q;
q_proposed = q;

Q_basis_proposed = cell(1,length(t));
Q_basis_proposed{1} = Q_proposed(:,:,1);

lambda_order = Polynomial.integer_grid([1,1,1]);
lambda_max_order = max(lambda_order,[],1);

xiQp_int = zeros(2,2);

x_region = 1 - (xp-q)'*inv(Q)*(xp-q);
x_region_order = x_region.order;
x_region_max_order = max(x_region_order,[],1);
E_region = Polynomial.expand_matrix(x_region_max_order, x_region_order);
x_region_coeff = E_region*x_region.coeff;

B_lambda = Polynomial.bernstein_transform_matrix(lambda_max_order);
E_lambda = Polynomial.expand_matrix(lambda_max_order, lambda_order);

P = Polynomial.product_matrix(x_region_max_order, lambda_max_order);
prod_max_order = x_region_max_order + lambda_max_order;
prod_full_order = Polynomial.integer_grid(prod_max_order);

xi_order = Polynomial.integer_grid([1,1,1]);
xi_max_order = [1,1,1];

const_max_order = max(xi_max_order, prod_max_order);
B_const = Polynomial.bernstein_transform_matrix(const_max_order);

E_xi = Polynomial.expand_matrix(const_max_order, xi_order);
E_affine = Polynomial.expand_matrix(const_max_order, [zeros(1,3); eye(3)]);
E_prod = Polynomial.expand_matrix(const_max_order, prod_full_order);

% decision variable: [b(1); b(2); a(nx); d(nw); lambda_1; lambda_2];
clear prob
prob.solver = 'linprog';
prob.options = optimoptions('linprog', 'Display', 'none');

prob.f = ...
    [1;...
    -1;...
    zeros(2,1);...
    0;...
    zeros(size(lambda_order,1),1)];

prob.const1_Amat = ...
    [-E_affine(:,1),...
    zeros(size(E_affine,1),1),...
    -E_affine(:,2:end),...
    E_prod*P*kron(x_region_coeff, E_lambda)];

prob.lambda1_Amat = ...
    [zeros(size(B_lambda,1), 2+2+1),...
    -E_lambda];

prob.const2_Amat = ...
    [zeros(size(E_affine,1),1),...
    E_affine(:,1),...
    E_affine(:,2:end),...
    E_prod*P*kron(x_region_coeff, E_lambda)];

n_const = size(B_const,1);
n_lambda = size(B_lambda,1);
prob.Aineq = zeros( n_const+n_const+n_lambda, 2+2+1+n_lambda );
prob.bineq = zeros( n_const+n_const+n_lambda, 1);

tic
for k = 1:length(t)-1
    % nominal values
    dt_ = t(k+1)-t(k);
    Q_ = Q_proposed(:,:,k);
    q_ = q_proposed(:,k);
    A0_ = dfdx(q_, 0);
    D0_ = dfdw(q_, 0);
    f0_ = f(q_, 0);
    xMax_ = sqrt(sum(Q_^(1/2).^2,2));
    
    % nonlinearity at the current instance
    xip_ = fp - (f0_ + A0_*(xp - q_) + D0_*(wp));
        
    % principal axes of nonlinearity
    xiQp_ = xip_ * xip_';
    for i = 1:2
        for j = 1:2
            xiQp_ij = xiQp_(i,j);
            xiQp_ij_int_order = xiQp_ij.order + 1;
            xiQp_ij_int_coeff = xiQp_ij.coeff./prod(xiQp_ij_int_order,2);
            xiQp_ij_int = sum(xiQp_ij_int_coeff.*prod(repmat([xMax_',wMax], [size(xiQp_ij_int_order,1),1]).^xiQp_ij_int_order,2))...
                - sum(xiQp_ij_int_coeff.*prod(repmat(-[xMax_',wMax], [size(xiQp_ij_int_order,1),1]).^xiQp_ij_int_order,2));
            xiQp_int(i,j) = xiQp_ij_int;
        end
    end
    [V_,~] = eig(xiQp_int);
    
    % compute nonlinerity bounder
    xi_prime_ = V_*xip_;
    
    Aprime_ = zeros(2,2);
    Dprime_ = zeros(2,1);
    bprime_ = zeros(2,2);
    
    scale_ = [xMax_; wMax];
    shift_ = [q_; 0];
    T_lambda_ = Polynomial.domain_transform_matrix(lambda_max_order, scale_, shift_);
    T_const_ = Polynomial.domain_transform_matrix(const_max_order, scale_, shift_);
    
    x_region_ = 1 - (xp-q_)'*inv(Q_)*(xp-q_);
    x_region_coeff_ = E_region*x_region_.coeff;
    
    prob.const1_Amat(:,5+(1:n_lambda)) = E_prod*P*kron(x_region_coeff_, E_lambda);
    prob.const2_Amat(:,5+(1:n_lambda)) = E_prod*P*kron(x_region_coeff_, E_lambda);

    prob.Aineq(1:n_const,:) = (B_const\T_const_)*prob.const1_Amat;
    prob.Aineq(n_const+(1:n_const),:) = (B_const\T_const_)*prob.const2_Amat;
    prob.Aineq(2*n_const+(1:n_lambda),:) = (B_lambda\T_lambda_)*prob.lambda1_Amat;
    
    for j = 1:2
        xi_coeff_ = xi_prime_(j).coeff;
        
        prob.bineq(1:n_const) = -(B_const\T_const_)*(E_xi*xi_coeff_);
        prob.bineq(n_const+(1:n_const)) = (B_const\T_const_)*(E_xi*xi_coeff_);
        
        opt_ = linprog(prob);
        Aprime_(j,:) = opt_(3:4);
        Dprime_(j,:) = opt_(5);
        bprime_(j,:) = 0.5*opt_(1:2);
    end
    
    % matrix manipulation
    A_ = A0_ + V_\Aprime_;
    D_ = [D0_ + V_\Dprime_, inv(V_)];
    c_ = sum(bprime_,2)/2;
    xiMax_ = abs(diff(bprime_,[],2)/2);
    wMax_ = [wMax; xiMax_];
    TF_ = expm(A_*dt_);
    
    % compute sets due to disturbances
    Nk_ = zeros(2,2,1+2);
    for j = 1:size(D_,2)
        M_ = (D_(:,j)*D_(:,j)');
        N_ = lyap(A_, expm(-A_*dt_)*M_*expm(-A_'*dt_) - M_);
        Nk_(:,:,j) = dt_*(wMax_(j)*wMax_(j))*N_;
    end
    
    % prepare for basis
    N_basis_prev = size(Q_basis_proposed{k},3);
    N_basis = N_basis_prev;
    N_basis = N_basis + 1 + 2;
    Q_basis_ = zeros(2,2,N_basis);
    Q_basis_(:,:,1:N_basis_prev) = Q_basis_proposed{k};
    
    % effect of disturbance
    Q_basis_(:,:,N_basis_prev+(1:(1+2))) = Nk_;
    for i = 1:size(Q_basis_,3)
        Q_basis_(:,:,i) = TF_ * Q_basis_(:,:,i) * TF_';
    end
    
    % composition
    Q_basis_proposed{k+1} = Q_basis_;
    Q_proposed(:,:,k+1) = Math.Ellipsoid.MVOE( Q_basis_ );
    q_proposed(:,k+1) = q_ + dt_*(A_*q_ + f0_ - A0_*q_ + V_\c_);
end
ctime_proposed = toc;
ctime_proposed
%% 2D visualization
nPts = 10;
nPts2 = 200;
figure(11)
set(gcf, 'position', [681,559,480,420])
cla; hold on; grid on; axis equal; axis tight;
for k = round(linspace(1,length(t),4))
    if k > 1
        %%% Stochastic propagation
        P_ = P_gmm(:,:,k);
        %%%%% draw 4.5-sigma line of sum of Gaussians
        tmp = sys.xN(:,k) + 4.5*P_^(1/2)*Math.Sphere(1,nPts).x;
        h1 = plot(tmp(1,:), tmp(2,:), 'x', 'color', [237,145,33]/255, 'linewidth', 1);
        tmp = sys.xN(:,k) + 4.5*P_^(1/2)*Math.Sphere(1,nPts2).x;
        plot(tmp(1,:), tmp(2,:), '-', 'color', [237,145,33]/255, 'linewidth', 1);
        
        %%% Linearization-based method
        x_ = sys.xN(:,k) + Q_lin(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h2 = plot(x_(1,:), x_(2,:), 'k.', 'linewidth', 1);
        x_ = sys.xN(:,k) + Q_lin(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), 'k-', 'linewidth', 1);
        
        %%% SOS-program-based method
        x_ = sys.xN(:,k) + Q_sos(:,:,k)^(1/2)*Math.Sphere(1,nPts).x;
        h3 = plot(x_(1,:), x_(2,:), '^', 'color', 'r', 'linewidth', 1);
        x_ = sys.xN(:,k) + Q_sos(:,:,k)^(1/2)*Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', 'r', 'linewidth', 1);
        
        %%% Nonlinear-optimization-based method
        x_ = sys.xN(:,k) + Q_nonlin(:,:,k)^(1/2)*Math.Sphere(1,nPts).x;
        h4 = plot(x_(1,:), x_(2,:), 's', 'color', 'b', 'linewidth', 1);
        x_ = sys.xN(:,k) + Q_nonlin(:,:,k)^(1/2)*Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', 'b', 'linewidth', 1);
        
        %%% Proposed method
        x_ = sys.xN(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h5 = plot(x_(1,:), x_(2,:), 'o', 'color', [65,169,76]/255, 'linewidth', 1);
        x_ = sys.xN(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', [65,169,76]/255, 'linewidth', 1);
    end
    
    %%% HJB eqn
    x_ = sys.xN(:,k) + X{k};
    h6 = plot(x_(1,:), x_(2,:), 'k-', 'linewidth', 2);
    
    text(sys.xN(1,k), sys.xN(2,k), ['$t=',num2str(t(k)),'$ s'],...
        'horizontalalignment', 'center', 'fontsize', 14)
end
legend([h6,h2,h3,h4,h1,h5],...
    'HJB PDE',...
    'Linearization',...
    'SOS Program',...
    'Nonlinear optimization',...
    'Uncertainty propagation',...
    'Proposed',...
    'location', 'southeast')
xlabel('$x_1$')
ylabel('$x_2$')
set(gca, 'xlim', [0.75, 1.25])
set(gca, 'ylim', [0.83, 1.3])

%% 3D visualization
F_proposed = cell(1,length(t));
F_sos = cell(1,length(t));
for k = 1:length(t)
    F_proposed{k} = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,500).x;
    F_sos{k} = res_sos(end).step2(:,:,k)^(1/2)*Math.Sphere(1,500).x;
end

figure;
cla; hold on; grid on; axis tight;
set(gcf, 'position', [681,559,480,420])
h7 = Utils.plot_set(X, t, length(t), 'k', 0.3);
h9 = Utils.plot_set(F_sos, t, length(t), 'r', 0.3);
h8 = Utils.plot_set(F_proposed, t, length(t), [65,169,76]/255, 0.3);

for k = round(linspace(1,length(t),5))
    if k > 1
        %%% SOS
        x_ = Q_sos(:,:,k)^(1/2)*Math.Sphere(1,500).x;
        h3 = plot3(x_(1,:), t(k)*ones(size(x_(1,:))), x_(2,:), '-', 'color', 'r', 'linewidth', 2);
        %%% Proposed
        x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,500).x;
        h5 = plot3(x_(1,:), t(k)*ones(size(x_(1,:))), x_(2,:), '-', 'color', [65,169,76]/255, 'linewidth', 2);
    end
    %%% HJB eqn
    x_ = X{k};
    h6 = plot3(x_(1,:), t(k)*ones(size(x_(1,:))), x_(2,:), 'k-', 'linewidth', 2);
end
xlabel('$\tilde{x}_1$')
ylabel('$t$ [s]')
zlabel('$\tilde{x}_2$')
view([122,11])
xlim([-0.1,0.1])
zlim([-0.1,0.1])
ylim([t(1), t(end)])
set(gca, 'xdir', 'reverse')
camlight('left')
legend([h7,h9,h8],...
    'HJB PDE',...
    'SOS program',...
    'Proposed',...
    'location', 'northwest')
camlight('right')

%% Computation times
disp('Computation time of the methods')
disp(['1) Level-set: ', num2str(ctime_levelset), ' seconds.'])
disp(['2) Linearization-based: ', num2str(ctime_linearization), ' seconds.'])
disp(['3) SOS-program: ', num2str(ctime_sos), ' seconds.'])
disp(['4) Nonlinear-optimization-based: ', num2str(ctime_nonlinopt), ' seconds.'])
disp(['5) Uncertainty propagation: ', num2str(ctime_stochastic), ' seconds.'])
disp(['6) Proposed: ', num2str(ctime_proposed), ' seconds.'])

%%% remove added path
rmpath(genpath('3rd_party/helperOC-master'))
rmpath(genpath('3rd_party/ToolboxLS'))
rmpath(genpath('3rd_party/SOSTOOLS'))
rmpath(genpath('3rd_party/SeDuMi_1_3'))
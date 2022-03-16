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
f = @(x,w,t)...
    [(w+3).*x(1,:).*((1-t)-x(2,:));...
    (w+3).*x(2,:).*(x(1,:)-(1-t))];
dfdx = @(x,w,t)...
    [(w+3)*((1-t)-x(2)), -(w+3)*x(1);...
    (w+3)*x(2), (w+3)*(x(1)-(1-t))];
dfdw = @(x,w,t)...
    [x(1)*((1-t)-x(2));...
    x(2)*(x(1)-(1-t))];

q = [1.2; 1.1];
Q = 0.05^2*eye(2);
wMax = 0.2;

% polynomial expression of dynamics
xp = [Polynomial(1,[1,0,0,0]); Polynomial(1,[0,1,0,0])];
wp = Polynomial(1,[0,0,1,0]);
tp = Polynomial(1,[0,0,0,1]);
fp = f(xp, wp, tp);
t = linspace(0,1,6);

% buffer
Q_proposed = zeros(2,2,length(t));
Q_proposed(:,:,1) = Q;
q_proposed = q;

Q_basis_proposed = cell(1,length(t));
Q_basis_proposed{1} = Q_proposed(:,:,1);

lambda_order = Polynomial.integer_grid([1,1,1,1]);
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

xi_order = Polynomial.integer_grid([1,1,1,1]);
xi_max_order = [1,1,1,1];

const_max_order = max(xi_max_order, prod_max_order);
B_const = Polynomial.bernstein_transform_matrix(const_max_order);


E_xi = Polynomial.expand_matrix(const_max_order, xi_order);
E_affine = Polynomial.expand_matrix(const_max_order, [zeros(1,4); [eye(3), zeros(3,1)]]);
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


%%
tic
for k = 1:length(t)-1

dt_ = t(k+1)-t(k);
Q_ = Q_proposed(:,:,k);
q_ = q_proposed(:,k);
A0_ = dfdx(q_, 0, t(k));
D0_ = dfdw(q_, 0, t(k));
f0_ = f(q_, 0, t(k));
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
        xiQp_ij_int = sum(xiQp_ij_int_coeff.*prod(repmat([xMax_',wMax,t(k+1)], [size(xiQp_ij_int_order,1),1]).^xiQp_ij_int_order,2))...
            - sum(xiQp_ij_int_coeff.*prod(repmat(-[xMax_',wMax,t(k)], [size(xiQp_ij_int_order,1),1]).^xiQp_ij_int_order,2));
        xiQp_int(i,j) = xiQp_ij_int;
    end
end
[V_,~] = eig(xiQp_int);

% compute nonlinerity bounder
xi_prime_ = V_*xip_;

Aprime_ = zeros(2,2);
Dprime_ = zeros(2,1);
bprime_ = zeros(2,2);

scale_ = [xMax_; wMax; 0.5*(t(k+1)-t(k))];
shift_ = [q_; 0; 0.5*(t(k+1) + t(k))];
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
    xi_coeff_ = Polynomial.expand_matrix(xi_max_order,xi_prime_(j).order)*xi_prime_(j).coeff;
    
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

%% 2D visualization
nPts = 10;
nPts2 = 200;
figure(11)
set(gcf, 'position', [681,559,480,420])
cla; hold on; grid on; axis equal; axis tight;
for k = round(linspace(1,length(t),11))
    if k > 1
%         %%% Stochastic propagation
%         P_ = P_gmm(:,:,k);
%         %%%%% draw 4.5-sigma line of sum of Gaussians
%         tmp = sys.xN(:,k) + 4.5*P_^(1/2)*Math.Sphere(1,nPts).x;
%         h1 = plot(tmp(1,:), tmp(2,:), 'x', 'color', [237,145,33]/255, 'linewidth', 1);
%         tmp = sys.xN(:,k) + 4.5*P_^(1/2)*Math.Sphere(1,nPts2).x;
%         plot(tmp(1,:), tmp(2,:), '-', 'color', [237,145,33]/255, 'linewidth', 1);
        
%         %%% Linearization-based method
%         x_ = sys.xN(:,k) + Q_lin(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         h2 = plot(x_(1,:), x_(2,:), 'k.', 'linewidth', 1);
%         x_ = sys.xN(:,k) + Q_lin(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         plot(x_(1,:), x_(2,:), 'k-', 'linewidth', 1);
        
%         %%% SOS-program-based method
%         x_ = sys.xN(:,k) + Q_sos(:,:,k)^(1/2)*Math.Sphere(1,nPts).x;
%         h3 = plot(x_(1,:), x_(2,:), '^', 'color', 'r', 'linewidth', 1);
%         x_ = sys.xN(:,k) + Q_sos(:,:,k)^(1/2)*Math.Sphere(1,nPts2).x;
%         plot(x_(1,:), x_(2,:), '-', 'color', 'r', 'linewidth', 1);
        
%         %%% Nonlinear-optimization-based method
%         x_ = sys.xN(:,k) + Q_nonlin(:,:,k)^(1/2)*Math.Sphere(1,nPts).x;
%         h4 = plot(x_(1,:), x_(2,:), 's', 'color', 'b', 'linewidth', 1);
%         x_ = sys.xN(:,k) + Q_nonlin(:,:,k)^(1/2)*Math.Sphere(1,nPts2).x;
%         plot(x_(1,:), x_(2,:), '-', 'color', 'b', 'linewidth', 1);
        
        %%% Proposed method
%         x_ = q_proposed(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h5 = plot(x_(1,:), x_(2,:), 'o', 'color', 'r', 'linewidth', 1);
%         x_ = q_proposed(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', 'r', 'linewidth', 1);
    end
    
%     %%% HJB eqn
%     x_ = sys.xN(:,k) + X{k};
%     h6 = plot(x_(1,:), x_(2,:), 'k-', 'linewidth', 2);
    
%     text(q_proposed(1,k), q_proposed(2,k), ['$t=',num2str(t(k)),'$ s'],...
%         'horizontalalignment', 'center', 'fontsize', 14)
end

%%
load('discrete.mat')
for k = round(linspace(1,length(t),11))
    if k > 1
%         %%% Stochastic propagation
%         P_ = P_gmm(:,:,k);
%         %%%%% draw 4.5-sigma line of sum of Gaussians
%         tmp = sys.xN(:,k) + 4.5*P_^(1/2)*Math.Sphere(1,nPts).x;
%         h1 = plot(tmp(1,:), tmp(2,:), 'x', 'color', [237,145,33]/255, 'linewidth', 1);
%         tmp = sys.xN(:,k) + 4.5*P_^(1/2)*Math.Sphere(1,nPts2).x;
%         plot(tmp(1,:), tmp(2,:), '-', 'color', [237,145,33]/255, 'linewidth', 1);
        
%         %%% Linearization-based method
%         x_ = sys.xN(:,k) + Q_lin(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         h2 = plot(x_(1,:), x_(2,:), 'k.', 'linewidth', 1);
%         x_ = sys.xN(:,k) + Q_lin(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         plot(x_(1,:), x_(2,:), 'k-', 'linewidth', 1);
        
%         %%% SOS-program-based method
%         x_ = sys.xN(:,k) + Q_sos(:,:,k)^(1/2)*Math.Sphere(1,nPts).x;
%         h3 = plot(x_(1,:), x_(2,:), '^', 'color', 'r', 'linewidth', 1);
%         x_ = sys.xN(:,k) + Q_sos(:,:,k)^(1/2)*Math.Sphere(1,nPts2).x;
%         plot(x_(1,:), x_(2,:), '-', 'color', 'r', 'linewidth', 1);
        
%         %%% Nonlinear-optimization-based method
%         x_ = sys.xN(:,k) + Q_nonlin(:,:,k)^(1/2)*Math.Sphere(1,nPts).x;
%         h4 = plot(x_(1,:), x_(2,:), 's', 'color', 'b', 'linewidth', 1);
%         x_ = sys.xN(:,k) + Q_nonlin(:,:,k)^(1/2)*Math.Sphere(1,nPts2).x;
%         plot(x_(1,:), x_(2,:), '-', 'color', 'b', 'linewidth', 1);
        
        %%% Proposed method
        x_ = q_proposed(:,k) + Q_discrete(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h5 = plot(x_(1,:), x_(2,:), 'o', 'color', 'b', 'linewidth', 1);
        x_ = q_proposed(:,k) + Q_discrete(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', 'b', 'linewidth', 1);
    end
    
%     %%% HJB eqn
%     x_ = sys.xN(:,k) + X{k};
%     h6 = plot(x_(1,:), x_(2,:), 'k-', 'linewidth', 2);
    
%     text(q_proposed(1,k), q_proposed(2,k), ['$t=',num2str(t(k)),'$ s'],...
%         'horizontalalignment', 'center', 'fontsize', 14)
end

% legend([h6,h2,h3,h4,h1,h5],...
%     'HJB PDE',...
%     'Linearization',...
%     'SOS Program',...
%     'Nonlinear optimization',...
%     'Uncertainty propagation',...
%     'Proposed',...
%     'location', 'southeast')
xlabel('$x_1$')
ylabel('$x_2$')
% set(gca, 'xlim', [0.75, 1.25])
% set(gca, 'ylim', [0.83, 1.3])

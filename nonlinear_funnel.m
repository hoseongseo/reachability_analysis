clear; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))

addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)


%% Dynamics / initial set / disturbance
% original equation
wMax = 0.1;
f = @(x,w) [(wMax*w+3).*x(1,:).*(1-x(2,:)); (wMax*w+3).*x(2,:).*(x(1,:)-1)];
% dfdx = @(x,w) [(wMax*w+3)*(1-x(2)), -(wMax*w+3)*x(1);...
%                (wMax*w+3)*x(2), (wMax*w+3)*(x(1)-1)];
% dfdw = @(x,w) wMax*[x(1)*(1-x(2)); x(2)*(x(1)-1)];

% transformed to local coordinate
g = @(xb,w,x0) f(xb+x0,w) - f(x0,0);
% dgdx = @(xb,w,x0) dfdx(xb+x0,w);
% dgdw = @(xb,w,x0) dfdw(xb+x0,w);

Q0 = diag([0.05; 0.05])^2; % initial set

%% Reference generation
t = linspace(0,1,101);

%%
sys0 = Dynamics.LotkaVolterraNominal;
sys = Dynamics.LotkaVolterra(wMax, sys0, [1.2; 1.1], t);

nGrid = [201, 201];
minGrid = [-0.06, -0.06];
maxGrid = [0.06, 0.06];
gr = createGrid(minGrid, maxGrid, nGrid);
V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) - 1;
X0 = getLevelSet(gr, V0, 0.0);

hjb_equation = HJBequation(sys, gr);
V = zeros([size(V0), length(t)]);
V(:,:,1) = V0;
for i = 1:length(t)-1
    V(:,:,i+1) = hjb_equation.solve(V(:,:,i), t(i), t(i+1));
end

X = cell(1,length(t));
X{1} = X0;
for i = 2:length(t)
    X{i} = getLevelSet(gr, V(:,:,i), 0.0);
end

%% Funnel (SOS Program)
tic
[Q_sos1, Q_sos2] = funnel_sos(...
    sys,... % system (polynomial dynamics)
    t,...   % time
    Q0,...  % initial shape matrix
    5.0,... % initial guess parameter
    21,...  % maximum iteration
    5e-4);  % convergence tolerance
toc

%% 
S = Utils.Sphere(1,200);
F_sos1 = zeros([size(S.x), length(t)]);
F_sos2 = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_sos1(:,:,i) = Q_sos1(:,:,i)^(1/2) * S.x;
    F_sos2(:,:,i) = Q_sos2(:,:,i)^(1/2) * S.x;
end

%% plot
figure;
cla; hold on; grid on;
plotSet(X, t, length(t), 'k', 0.7);
plotSet(F_sos1, t, length(t), 'g', 0.5);
% plotSet(F_sos2, t, length(t), 'b', 0.5);
% plotSet(xChar(:,:,local_idx), t(local_idx), length(local_idx), 'k', 0.5);
view([128,11])
camlight left
camlight right
xlabel('$x_1$')
ylabel('$t$ [s]')
zlabel('$x_2$')
ax = gca;
ax.XLabel.Position = [-0.0056, 1.0799, -0.0657];
ax.YLabel.Position = [0.0625, 0.4441, -0.0690];
% legend([h1,h2], '$\underline{\mathcal{X}}(t)$', '$\mathcal{E}(Q_x(t))$', 'location', 'northwest')

%% Proposed
tk = t(local_idx);
Ak = ANominal(:,:,local_idx);
Dk = DNominal(:,:,local_idx);
xk = x0(:,local_idx);
wMag = (max(W) - min(W))*0.5;

tic
Q_frs_lin = ellipsoidal_reachability(tk, Ak, Dk, xk, wMag, Q0, 'determinant', false, g);
toc
%%
Q_frs_nonlin = ellipsoidal_reachability(tk, Ak, Dk, xk, wMag, Q0, 'determinant', true, g);

%%
tic
Q_frs_nonlin2 = ellipsoidal_reachability2(tk, Ak, Dk, xk, wMag, Q0, 'determinant', true, g);
toc

%% DIRTREL
Nk = length(tk);

W = wMax^2; % disturbance bound
Ek = zeros(2,2,Nk);
Ek(:,:,1) = Q0;
Hk = zeros(2,1,Nk);
Hk(:,:,1) = 1e-20*ones(2,1);

E_basis = cell(1,Nk);
E_basis{1} = Q0;

for i = 1:Nk-1 % discrete version
    dt_ = tk(i+1) - tk(i);
%     E_ = Ek(:,:,i);
    H_ = Hk(:,:,i);
    A_ = eye(2) + dt_*Ak(:,:,i);
    D_ = dt_*Dk(:,:,i);
    
%     Enew_ = A_*E_*A_' + A_*H_*D_' + D_*H_'*A_' + D_*W*D_';
%     Enew_ = Math.Ellipsoid.MVOE( cat(3, A_*E_*A_', A_*H_*D_', D_*H_'*A_', D_*W*D_') );
    E_basis_prev_ = E_basis{i};
%     E_ = Math.Ellipsoid.MVOE( E_basis_prev_ );
    
    E_basis_new_ = zeros(2,2,size(E_basis_prev_,3)+3);
    E_basis_new_(:,:,1:size(E_basis_prev_,3)) = E_basis_prev_;
    for k = 1:size(E_basis_prev_,3)
        E_basis_new_(:,:,k) = A_ * E_basis_new_(:,:,k) * A_';
    end
    E_basis_new_(:,:,end-2) = A_*H_*D_';
    E_basis_new_(:,:,end-1) = D_*H_'*A_';
    E_basis_new_(:,:,end) = D_*W*D_';
    
    Hnew_ = A_*H_ + D_*W;
%     Enew_ = Math.Ellipsoid.MVOE(E_basis_new_);
    Enew_ = sum(E_basis_new_,3);
    
    Ek(:,:,i+1) = Enew_;
    Hk(:,:,i+1) = Hnew_;
    E_basis{i+1} = E_basis_new_;
end

tic
Q_frs_sum = ellipsoidal_reachability(tk, Ak, Dk, xk, wMag, Q0, 'sum', false, g);
toc


%%
figure(21)
cla; hold on; grid on;
% for i = 1:size(x0,2)
%     h5 = plot( xNominal(1,:)' + squeeze(xChar(1,i,:)), xNominal(2,:)' + squeeze(xChar(2,i,:)),...
%         'color', [0.5,0.5,0.5], 'linewidth', 0.1);
% end
for i = 1:size(x1,2)
    h5 = plot( x0(1,:)' + squeeze(xChar1(1,i,:)), x0(2,:)' + squeeze(xChar1(2,i,:)),...
        'color', [0.5,0.5,0.5], 'linewidth', 0.1);
end
% plot(xNominal(1,:), xNominal(2,:), 'k')
for k = round(linspace(1,length(local_idx),5))
    i = local_idx(k);
    h7 = plot(x0(1,i) + xChar1(1,:,i), x0(2,i) + xChar1(2,:,i), 'm', 'linewidth', 2);
    h1 = plot(x0(1,i) + xChar(1,:,i), x0(2,i) + xChar(2,:,i), 'k', 'linewidth', 2);
    
%     tmp = xNominal(:,i) + Q_frs_lin(:,:,k)^(1/2)*Math.Sphere(1,200).x;
%     h2 = plot(tmp(1,:), tmp(2,:), 'color', [0,128,0]/256, 'linewidth', 2);
    
%     tmp = xNominal(:,i) + Q_frs_nonlin(:,:,k)^(1/2)*Math.Sphere(1,200).x;
%     plot(tmp(1,:), tmp(2,:), 'b', 'linewidth', 2)
    
%     tmp = xNominal(:,i) + Q_frs_nonlin2(:,:,k)^(1/2)*Math.Sphere(1,200).x;
%     h3 = plot(tmp(1,:), tmp(2,:), 'b', 'linewidth', 2);
    
%     tmp = xNominal(:,i) + Q_funnel(:,:,k)^(1/2)*Math.Sphere(1,200).x;
%     h4 = plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2);
%     
%     tmp = xNominal(:,i) + Ek(:,:,k)^(1/2)*Math.Sphere(1,200).x;
%     h6 = plot(tmp(1,:), tmp(2,:), 'b-', 'linewidth', 2);
    
    text(x0(1,i), x0(2,i), ['$t=',num2str(t(i)),'$ s'],...
        'horizontalalignment', 'center', 'fontsize', 14)
end
axis tight;
set(gca, 'fontsize', 14)
xlabel('$x_1$')
ylabel('$x_2$')
% legend([h1,h2,h3,h4],...
%     'True set',...
%     'Proposed without $\mathcal{E}(Q_\Omega)$',...
%     'Proposed with $\mathcal{E}(Q_\Omega)$',...
%     'Funnel',...
%     'location', 'southeast')
% legend([h1,h2,h4],...
%     'True set',...
%     'Proposed',...
%     'Funnel',...
%     'location', 'southeast')
% legend([h1,h3,h4],...
%     'Reachable sets',...
%     'Proposed approximation',...
%     'Funnel',...
%     'location', 'southeast')
legend([h5,h7,h1],...
    'Solutions of characteristic eqns.',...
    'Collection of states at the same instance',...
    'Reachable sets',...
    'location', 'southeast')
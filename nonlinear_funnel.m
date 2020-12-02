clear; clc;

addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver

%% Dynamics / initial set / disturbance
% original equation
f = @(x,w) [(0.1*w+3).*x(1,:).*(1-x(2,:)); (0.1*w+3).*x(2,:).*(x(1,:)-1)];
dfdx = @(x,w) [(0.1*w+3)*(1-x(2)), -(0.1*w+3)*x(1); (0.1*w+3)*x(2), (0.1*w+3)*(x(1)-1)];
dfdw = @(x,w) 0.1*[x(1)*(1-x(2)); x(2)*(x(1)-1)];

% transformed to local coordinate
g = @(xb,w,x0) f(xb+x0,w) - f(x0,0); 
dgdx = @(xb,w,x0) dfdx(xb+x0,w);
dgdw = @(xb,w,x0) dfdw(xb+x0,w);

Q = diag([0.05; 0.05])^2; % initial set
wgr = linspace(-1, 1, 11);

%% Reference generation
t = linspace(0,1,10001); % this should be as precise as possible
dt = t(2) - t(1);
N = length(t);

xInit = [1.2; 1.1];

xNominal = zeros(2, length(t));
xNominal(:,1) = xInit;
for i = 1:N-1
    x_ = xNominal(:,i);
    f_ = f(x_, 0);
    xNominal(:,i+1) = x_ + dt*f_;
end

% linearized dynamics
ANominal = zeros(2,2,N);
DNominal = zeros(2,1,N);
for i = 1:N
    ANominal(:,:,i) = dfdx( xNominal(:,i), 0 );
    DNominal(:,:,i) = dfdw( xNominal(:,i), 0 );
end

local_idx = 1:100:10001; % time stamp

%% Characteristic equation
x0 = Q^(1/2)*Math.Sphere(1,50).x;
p0 = 2*Q\x0;
z0 = zeros(1,size(x0,2));
xChar = zeros(length(xInit),size(x0,2),length(t));
pChar = zeros(length(xInit),size(p0,2),length(t));
zChar = zeros(1,size(x0,2),length(t));
xChar(:,:,1) = x0;
pChar(:,:,1) = p0;
zChar(:,:,1) = z0;
for k = 1:size(x0,2)
    for i = 1:length(t)-1   
        xn_ = xNominal(:,i);
        x_ = xChar(:,k,i);
        p_ = pChar(:,k,i);
        z_ = zChar(:,k,i);
        
        HAll_ = zeros(1,length(wgr));
        for j = 1:length(wgr)
            w_ = wgr(j);
            f_ = g(x_, w_, xn_);
            HAll_(j) = dot(p_, f_);
        end
        [H_, ind_] = max(HAll_);
        wOpt_ = wgr(ind_);
        
        Hp_ = g(x_,wOpt_,xn_);
        Hx_ = dgdx(x_,wOpt_,xn_)'*p_;
        
        dx_ = Hp_;
        dp_ = -Hx_;
        dz_ = dot(p_, Hp_) - H_;
        xChar(:,k,i+1) = x_ + dt*dx_;
        pChar(:,k,i+1) = p_ + dt*dp_;
        zChar(:,k,i+1) = z_ + dt*dz_;
    end
end

x1 = Q^(1/2)*Math.Sphere(1,8).x;
p1 = 2*Q\x1;
z1 = zeros(1,size(x1,2));
xChar1 = zeros(length(xInit),size(x1,2),length(t));
pChar1 = zeros(length(xInit),size(p1,2),length(t));
zChar1 = zeros(1,size(x1,2),length(t));
xChar1(:,:,1) = x1;
pChar1(:,:,1) = p1;
zChar1(:,:,1) = z1;
for k = 1:size(x1,2)
    for i = 1:length(t)-1   
        xn_ = xNominal(:,i);
        x_ = xChar1(:,k,i);
        p_ = pChar1(:,k,i);
        z_ = zChar1(:,k,i);
        
        HAll_ = zeros(1,length(wgr));
        for j = 1:length(wgr)
            w_ = wgr(j);
            f_ = g(x_, w_, xn_);
            HAll_(j) = dot(p_, f_);
        end
        [H_, ind_] = max(HAll_);
        wOpt_ = wgr(ind_);
        
        Hp_ = g(x_,wOpt_,xn_);
        Hx_ = dgdx(x_,wOpt_,xn_)'*p_;
        
        dx_ = Hp_;
        dp_ = -Hx_;
        dz_ = dot(p_, Hp_) - H_;
        xChar1(:,k,i+1) = x_ + dt*dx_;
        pChar1(:,k,i+1) = p_ + dt*dp_;
        zChar1(:,k,i+1) = z_ + dt*dz_;
    end
end

%% Funnel
tic
Q_sos = funnel_sos( t(local_idx), xNominal(:,local_idx), g, Q, 1, 5.0, 21, 5e-4, xChar(:,:,local_idx) );
toc

%% plot
S = Utils.Sphere(1,200);
F_sos = zeros([size(S.x), length(local_idx)]);
for i = 1:length(local_idx)
    F_sos(:,:,i) = Q_sos(:,:,i)^(1/2) * S.x;
end

figure;
cla; hold on; grid on;
% h1 = plotSet(X, t, length(t), 'k', 0.7);
plotSet(F_sos, t(local_idx), length(local_idx), 'g', 0.5);
plotSet(xChar(:,:,local_idx), t(local_idx), length(local_idx), 'k', 0.5);
view([128,11])
camlight left
camlight right
xlabel('$x_1$')
ylabel('$t$ [s]')
zlabel('$x_2$')
% ax = gca;
% ax.XLabel.Position = [-0.1725, 2.1516, -3.3917];
% ax.YLabel.Position = [2.2098, 0.9221, -3.5364];
% legend([h1,h2], '$\underline{\mathcal{X}}(t)$', '$\mathcal{E}(Q_x(t))$', 'location', 'northwest')
%% Proposed
tk = t(local_idx);
Ak = ANominal(:,:,local_idx);
Dk = DNominal(:,:,local_idx);
xk = xNominal(:,local_idx);
wMag = (max(W) - min(W))*0.5;

tic
Q_frs_lin = ellipsoidal_reachability(tk, Ak, Dk, xk, wMag, Q, 'determinant', false, g);
toc
%%
Q_frs_nonlin = ellipsoidal_reachability(tk, Ak, Dk, xk, wMag, Q, 'determinant', true, g);

%%
tic
Q_frs_nonlin2 = ellipsoidal_reachability2(tk, Ak, Dk, xk, wMag, Q, 'determinant', true, g);
toc

%% DIRTREL
Nk = length(tk);

W = wMax^2; % disturbance bound
Ek = zeros(2,2,Nk);
Ek(:,:,1) = Q;
Hk = zeros(2,1,Nk);
Hk(:,:,1) = 1e-20*ones(2,1);

E_basis = cell(1,Nk);
E_basis{1} = Q;

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
Q_frs_sum = ellipsoidal_reachability(tk, Ak, Dk, xk, wMag, Q, 'sum', false, g);
toc


%%
figure(21)
cla; hold on; grid on;
% for i = 1:size(x0,2)
%     h5 = plot( xNominal(1,:)' + squeeze(xChar(1,i,:)), xNominal(2,:)' + squeeze(xChar(2,i,:)),...
%         'color', [0.5,0.5,0.5], 'linewidth', 0.1);
% end
for i = 1:size(x1,2)
    h5 = plot( xNominal(1,:)' + squeeze(xChar1(1,i,:)), xNominal(2,:)' + squeeze(xChar1(2,i,:)),...
        'color', [0.5,0.5,0.5], 'linewidth', 0.1);
end
% plot(xNominal(1,:), xNominal(2,:), 'k')
for k = round(linspace(1,length(local_idx),5))
    i = local_idx(k);
    h7 = plot(xNominal(1,i) + xChar1(1,:,i), xNominal(2,i) + xChar1(2,:,i), 'm', 'linewidth', 2);
    h1 = plot(xNominal(1,i) + xChar(1,:,i), xNominal(2,i) + xChar(2,:,i), 'k', 'linewidth', 2);
    
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
    
    text(xNominal(1,i), xNominal(2,i), ['$t=',num2str(t(i)),'$ s'],...
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
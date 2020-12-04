clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

% given initial set
Q0 = diag([0.05; 0.05])^2;

% disturbance bound
wMax = 0.1;

% sphere
S = Utils.Sphere(1,200);

%%% system (polynomial dynamics)
t = linspace(0,1,101);
sys0 = Dynamics.LotkaVolterraNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
sys = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin


%%% Funnel (proposed, LTV)
Q_ltv = funnel_ltv(sys, t, Q0);
F_ltv = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_ltv(:,:,i) = Q_ltv(:,:,i)^(1/2) * S.x;
end

%%% Funnel (proposed, nonlinear)
[res_lp, cost_lp, rate_lp] = funnel_nonlinear_lp(sys, t, Q0, [3,3]);
F_lp = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_lp(:,:,i) = res_lp(end).Q(:,:,i)^(1/2) * S.x;
end

%%% Funnel (SOS Program)
[res_sos, cost_sos, rate_sos] = funnel_nonlinear_sos(sys, t, Q0);
F_sos = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_sos(:,:,i) = res_sos(end).step2(:,:,i)^(1/2) * S.x;
end

%%% HJB equation
% grid for HJB equation
nGrid = [201, 201];
minGrid = [-0.06, -0.06];
maxGrid = [0.06, 0.06];
gr = createGrid(minGrid, maxGrid, nGrid);
V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) - 1;
X0 = Utils.get_level_set(gr, V0, 0.0);

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


%%%%% comparison
figure;
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.5);
h2 = Utils.plot_set(F_lp, t, length(t), 'r', 0.3);
for i = round(linspace(1,length(t),7))
    tmp = F_sos(:,:,i);
    h3 = plot3(tmp(1,:), t(i)*ones(1,size(tmp,2)), tmp(2,:), 'color', [76,187,23]/255, 'linewidth', 2);
    tmp = F_lp(:,:,i);
    plot3(tmp(1,:), t(i)*ones(1,size(tmp,2)), tmp(2,:), 'r--', 'linewidth', 2)
end
view([128,11])
camlight left
camlight right
xlabel('$x_1$')
ylabel('$t$ [s]')
zlabel('$x_2$')
ax = gca;
ax.XLabel.Position = [-0.0056, 1.0799, -0.0657];
ax.YLabel.Position = [0.0625, 0.4441, -0.0690];
legend([h1,h2,h3],...
    '$\mathcal{X}(t)$',....
    '$\mathcal{F}(t)$ (Proposed)',...
    '$\mathcal{F}(t)$ (SOS program)',...
    'location', 'northeast')
title('$\textbf{Comparision of funnels and reachable set}$')

figure;
cla; hold on; grid on;
for i = round(linspace(1,length(t),5))
    tmp = sys.xN(:,i) + F_sos(:,:,i);
    h3 = plot(tmp(1,:), tmp(2,:), 'color', [76,187,23]/255, 'linewidth', 2);
    
    tmp = sys.xN(:,i) + F_lp(:,:,i);
    h2 = plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2);
    
    tmp = sys.xN(:,i) + X{i};
    h1 = plot(tmp(1,:), tmp(2,:), 'k', 'linewidth', 2);
    
    text(sys.xN(1,i), sys.xN(2,i), ['$t = ',num2str(t(i)),'$'],...
        'horizontalalignment', 'center', 'fontsize', 14)
end
xlabel('$x_1$')
ylabel('$x_2$')
legend([h1,h2,h3],...
    '$\mathcal{X}(t)$',....
    '$\mathcal{F}(t)$ (Proposed)',...
    '$\mathcal{F}(t)$ (SOS program)',...
    'location', 'southeast')
title('$\textbf{Funnels along the nominal trajectory}$')

figure;
subplot(2,1,1)
title('$\textbf{Convergence characteristics}$')
cla; hold on; grid on;
h1 = plot(cost_sos,'*-','color', [76,187,23]/255, 'linewidth', 2);
h2 = plot(cost_lp,'r*-', 'linewidth', 2);
legend([h1,h2],'Proposed', 'SOS program')
axis tight;
ylabel('$cost$')

subplot(2,1,2)
cla; hold on; grid on;
plot(rate_sos,'*-','color', [76,187,23]/255, 'linewidth', 2);
plot(rate_lp,'r*-', 'linewidth', 2)
axis tight;
ylabel('$rate$')
xlabel('Iterations')

%%% remove added path
rmpath(genpath('3rd_party/helperOC-master'))
rmpath(genpath('3rd_party/ToolboxLS'))
rmpath(genpath('3rd_party/SOSTOOLS'))
rmpath(genpath('3rd_party/SeDuMi_1_3'))
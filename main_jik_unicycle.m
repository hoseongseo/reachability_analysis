clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

set(groot, 'DefaultFigureColor', 'w')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
set(groot, 'DefaultAxesFontSize', 16)

% given initial set
Q0 = diag([0.1; 0.1; 0.1])^2;

% disturbance bound
% wMax = 0.1;

% sphere
S = Utils.Sphere(2,500);

%%% system (polynomial dynamics)
t = 0:0.05:2;
sys = Dynamics.JIK_unicycle(1.0, 0.5, 0.1);


%% Funnel (proposed, nonlinear)
% [res_lp, cost_lp, rate_lp] = funnel_nonlinear_lp(sys, t, Q0, [3,3]);
args = struct('max_iter', 30, 'ftol', 5e-4, 'plot_cost', true);

% [res_lp, cost_lp, rate_lp] = funnel_nonlinear_lp3(sys, t, Q0, [2,2], [2,2,2], [2,2], args);
tic
[res_lp, cost_lp, rate_lp] = funnel_nonlinear_lp4(sys, t, Q0, [2,2,2], [0,0,0,0], [2,2,2], args);
toc
F_lp = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_lp(:,:,i) = res_lp(end).Q(:,:,i)^(1/2) * S.x;
end


% %%
% %%% Funnel (SOS Program)
% args.max_iter = 10;
% args.rho = 3.0;
% args.ftol = 5e-4;
% args.plot_cost = true;
% tic;
% [res_sos, cost_sos, rate_sos] = funnel_nonlinear_sos(sys_poly, t, Q0, args);
% toc;
% F_sos = zeros([size(S.x), length(t)]);
% for i = 1:length(t)
%     F_sos(:,:,i) = res_sos(end).step2(:,:,i)^(1/2) * S.x;
% end

%%
%%% HJB equation
% grid for HJB equation
nGrid = [101, 101, 101];
minGrid = [-0.3, -0.6, -0.4];
maxGrid = [0.3, 0.6, 0.4];
gr = createGrid(minGrid, maxGrid, nGrid);
V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) + gr.xs{3}.*gr.xs{3}/Q0(3,3) - 1;
X0 = Utils.get_level_set(gr, V0, 0.0);

%%
% solve
hjb_equation = HJBequation(sys, gr);
V = zeros([size(V0), length(t)]);
V(:,:,:,1) = V0;
X = cell(1,length(t));
X{1} = X0;
figure(1)
cla; hold on; grid on; axis equal; axis tight;
xlim([minGrid(1), maxGrid(1)])
ylim([minGrid(2), maxGrid(2)])
zlim([minGrid(3), maxGrid(3)])
for i = 1:length(t)-1
    V(:,:,:,i+1) = hjb_equation.solve(V(:,:,:,i), t(i), t(i+1));
    X{i+1} = Utils.get_level_set(gr, V(:,:,:,i+1), 0.0);
    cla; hold on;
    plot3(X{i+1}(1,:), X{i+1}(2,:), X{i+1}(3,:), 'r.')
    plot3(F_lp(1,:,i+1), F_lp(2,:,i+1), F_lp(3,:,i+1), 'b.')
    title(['$t=',num2str(t(i+1)),'$'])
    xlabel('$\Delta x$ [m]')
    ylabel('$\Delta y$ [m]')
    zlabel('$\Delta \theta$ [rad]')
    view([-37,21])
    drawnow
end

%%

%%%%% comparison
figure;
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.5);
h2 = Utils.plot_set(F_lp, t, length(t), 'r', 0.3);
% h3 = Utils.plot_set(F_sos, t, length(t), 'r', 0.3);
% for i = round(linspace(1,length(t),7))
% %     tmp = F_sos(:,:,i);
% %     h3 = plot3(tmp(1,:), t(i)*ones(1,size(tmp,2)), tmp(2,:), 'color', [76,187,23]/255, 'linewidth', 2);
%     tmp = F_lp(:,:,i);
%     plot3(tmp(1,:), t(i)*ones(1,size(tmp,2)), tmp(2,:), 'r--', 'linewidth', 2)
% end
view([128,11])
camlight left
camlight right
xlabel('$x_1$')
ylabel('$t$ [s]')
zlabel('$x_2$')
ax = gca;
% ax.XLabel.Position = [-0.0056, 1.0799, -0.0657];
% ax.YLabel.Position = [0.0625, 0.4441, -0.0690];
% legend([h1,h2,h3],...
%     '$\mathcal{X}(t)$',....
%     '$\mathcal{F}(t)$ (Proposed)',...
%     '$\mathcal{F}(t)$ (SOS program)',...
%     'location', 'northeast')
title('$\textbf{Comparision of funnels and reachable set}$')

%%
figure;
cla; hold on; grid on;
for i = round(linspace(1,length(t),11))    
    tmp = F_lp(:,:,i);
    plot3(tmp(1,:), t(i)*ones(size(tmp(1,:))), tmp(2,:), 'r.', 'linewidth', 2)
end

%%
figure;
for k = 1:length(t)
subplot(3,1,1)
hold on; grid on; axis tight;
plot(t(k)*ones(size(F_lp(1,:,k))), F_lp(1,:,k), 'b.')
subplot(3,1,2)
hold on; grid on; axis tight;
plot(t(k)*ones(size(F_lp(2,:,k))), F_lp(2,:,k), 'b.')
subplot(3,1,3)
hold on; grid on; axis tight;
plot(t(k)*ones(size(F_lp(3,:,k))), F_lp(3,:,k), 'b.')
end
subplot(3,1,1)
ylabel('$\Delta x$ [m]')
subplot(3,1,2)
ylabel('$\Delta y$ [m]')
subplot(3,1,3)
ylabel('$\Delta \theta$ [rad]')
xlabel('$t$ [s]')

%%
figure;
cla; hold on; grid on;
for i = round(linspace(1,length(t),5))
    tmp = sys.xN(:,i) + F_sos(:,:,i);
    h3 = plot(tmp(1,:), tmp(2,:), 'color', [76,187,23]/255, 'linewidth', 2);
    
%     tmp = sys.xN(:,i) + F_lp(:,:,i);
%     h2 = plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2);
    
%     tmp = sys.xN(:,i) + X{i};
%     h1 = plot(tmp(1,:), tmp(2,:), 'k', 'linewidth', 2);
    
    text(sys.xN(1,i), sys.xN(2,i), ['$t = ',num2str(t(i)),'$'],...
        'horizontalalignment', 'center', 'fontsize', 14)
end
xlabel('$x_1$')
ylabel('$x_2$')
% legend([h1,h2,h3],...
%     '$\mathcal{X}(t)$',....
%     '$\mathcal{F}(t)$ (Proposed)',...
%     '$\mathcal{F}(t)$ (SOS program)',...
%     'location', 'southeast')
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
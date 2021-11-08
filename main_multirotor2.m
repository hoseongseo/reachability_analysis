clear; clc;

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
Q0 = (diag([0.03,0.03,0.03,0.03,0.03,0.03,0.01,0.01])^2);

% disturbance bound
wMax = [0.01; 0.01; 0.01];

% sphere
S = Utils.Sphere(2,500);

%%% system (polynomial dynamics)
t = linspace(0,1,51);
% t = linspace(0,0.1,11);
% rho = 1;
% u1 = rho*sin(2*pi*t) + 9.8;
% u2 = 0.1*sin(2*pi*t);
% u3 = 0.1*sin(2*pi*t);

xway = [0, 0.8];
yway = [0, 0.25, 0.4]; g = 9.8;
% yway = 0.0*[0, 0.25, 0.4]; g = 9.8;
% yway = -0.5*[0, 0.25, 0.4]; g = 9.8;
zway = [0, 0.1, 0.3];

x_IC = struct('vel', 0.0, 'acc', 0.0, 'jerk', 0.0, 'snap', []);
y_IC = struct('vel', 0.0, 'acc', 0.0, 'jerk', 0.0, 'snap', []);
z_IC = struct('vel', 0.0, 'acc', 0.0, 'jerk', 0.0, 'snap', []);
x_TC = struct('vel', [], 'acc', [], 'jerk', 0.0, 'snap', []);
y_TC = struct('vel', [], 'acc', [], 'jerk', 0.0, 'snap', []);
z_TC = struct('vel', [], 'acc', [], 'jerk', 0.0, 'snap', []);

x_tway = t(end);
y_tway = [0.55*t(end), 0.45*t(end)];
z_tway = [0.4*t(end), 0.6*t(end)];

traj_x = Planner.TrajectoryGeneration(xway, x_tway, 9, x_IC, x_TC);
traj_y = Planner.TrajectoryGeneration(yway, y_tway, 9, y_IC, y_TC);
traj_z = Planner.TrajectoryGeneration(zway, z_tway, 9, z_IC, z_TC);
ref_x = traj_x.at(t);
ref_y = traj_y.at(t);
ref_z = traj_z.at(t);

a = [ref_x(3,:)', ref_y(3,:)', ref_z(3,:)']';
j = [ref_x(4,:)', ref_y(4,:)', ref_z(4,:)']';
s = [ref_x(5,:)', ref_y(5,:)', ref_z(5,:)']';

[~, w, ~] = multirotor_flatness(a, j, s);
f = sqrt(sum((a + [0; 0; g]).^2,1));
% f = 9.8*ones(size(t));
% w = zeros(2,length(t));

sys0 = Dynamics.MultirotorNominal([0;0;0;0.0;0.0;0.0;0;0], [f; w(1,:); w(2,:)], t); % nominal dynamics with the initial condition
sys = Dynamics.Multirotor(sys0, wMax); % system shifted to the origin

sys0_x = Dynamics.MultirotorXNominal([0;0;0;0], [f; w(1,:); w(2,:)], t); % nominal dynamics with the initial condition
sys_x = Dynamics.MultirotorX(sys0_x, wMax(1));

sys0_y = Dynamics.MultirotorYNominal([0;0;0;0], [f; w(1,:); w(2,:)], t); % nominal dynamics with the initial condition
sys_y = Dynamics.MultirotorY(sys0_y, wMax(2));

sys0_z = Dynamics.MultirotorZNominal([0;0;0;0], [f; w(1,:); w(2,:)], t); % nominal dynamics with the initial condition
sys_z = Dynamics.MultirotorZ(sys0_z, wMax(3));

args.uN = [f; w(1,:); w(2,:)];
args.internal_composition = false;
args.composition_method = "MVOE";
args.eps = 1e-10;
Q_ltv = funnel_ltv(sys, t, Q0, args);

Q_ltv_x = funnel_ltv(sys_x, t, Q0([1,4,7,8],[1,4,7,8]), args);
Q_ltv_y = funnel_ltv(sys_y, t, Q0([2,5,7,8],[2,5,7,8]), args);
Q_ltv_z = funnel_ltv(sys_z, t, Q0([3,6,7,8],[3,6,7,8]), args);

figure(1)
hold on; grid on; axis equal;
plot3(sys0.xN(1,:), sys0.xN(2,:), sys0.xN(3,:), 'b', 'linewidth', 2)
plot3(ref_x(1,:), ref_y(1,:), ref_z(1,:), 'r--', 'linewidth', 2)
for k = round(linspace(1,length(t),11))
    Qtmp = Q_ltv(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + Qtmp(1:3,1:3) * Math.Sphere(2,300).x;
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'b.')
    
    Qtmp = Q_ltv_x(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [Qtmp(1)*[-1,1]; zeros(2,2)];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'g*-', 'linewidth', 2)
    
    Qtmp = Q_ltv_y(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [zeros(1,2); Qtmp(1)*[-1,1]; zeros(1,2)];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'g*-', 'linewidth', 2)
        
    Qtmp = Q_ltv_z(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [zeros(2,2); Qtmp(1)*[-1,1]];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'g*-', 'linewidth', 2)
    
    drawnow
end

view([28,28])
% view([-90,90])
% view([0,0])
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
zlabel('$x_3$ [m]')

%%
%%% Funnel (proposed, nonlinear)
% [res_lp, cost_lp, rate_lp] = funnel_nonlinear_lp(sys, t, Q0, [3,3]);
args = struct('max_iter', 1, 'ftol', 5e-4, 'plot_cost', false,...
    'grid', []);
args.internal_composition = false;
args.composition_method = "MVOE";
args.eps = 1e-10;
args.uN = [f; w(1,:); w(2,:)];
% args.uN = zeros(3,length(t));

%%
tic
res_lp_x = funnel_nonlinear_lp4(sys_x, t, Q0([1,4,7,8],[1,4,7,8]), [2,2,2,2], zeros(1,5), [2,2,2,2], args);
ctime_all_lp = toc;

%%
tic
res_lp_y = funnel_nonlinear_lp4(sys_y, t, Q0([2,5,7,8],[2,5,7,8]), [2,2,2,2], zeros(1,5), [2,2,2,2], args);
ctime_all_lp = toc;

%%
tic
res_lp_z = funnel_nonlinear_lp4(sys_z, t, Q0([3,6,7,8],[3,6,7,8]), [2,2,2,2], zeros(1,5), [2,2,2,2], args);
ctime_all_lp = toc;

%%
figure(1)
hold on; grid on; axis equal;
plot3(sys0.xN(1,:), sys0.xN(2,:), sys0.xN(3,:), 'b', 'linewidth', 2)
% plot3(ref_x(1,:), ref_y(1,:), ref_z(1,:), 'r--', 'linewidth', 2)
for k = round(linspace(1,length(t),11))
    Qtmp = Q_ltv(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + Qtmp(1:3,1:3) * Math.Sphere(2,300).x;
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'b.')
    
    Qtmp = Q_ltv_x(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [Qtmp(1)*[-1,1]; zeros(2,2)];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'g*-', 'linewidth', 2)
    
    Qtmp = res_lp_x(end).Q(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [Qtmp(1)*[-1,1]; zeros(2,2)];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'ro--', 'linewidth', 2)
    
    
    Qtmp = Q_ltv_y(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [zeros(1,2); Qtmp(1)*[-1,1]; zeros(1,2)];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'g*-', 'linewidth', 2)
    
    Qtmp = res_lp_y(end).Q(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [zeros(1,2); Qtmp(1)*[-1,1]; zeros(1,2)];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'ro--', 'linewidth', 2)
    
    Qtmp = Q_ltv_z(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [zeros(2,2); Qtmp(1)*[-1,1]];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'g*-', 'linewidth', 2)
    
    Qtmp = res_lp_z(end).Q(:,:,k)^(1/2);
    tmp = sys.xN(1:3,k) + [zeros(2,2); Qtmp(1)*[-1,1]];
    plot3(tmp(1,:), tmp(2,:), tmp(3,:), 'ro--', 'linewidth', 2)
    drawnow
end

view([28,28])
% view([-90,90])
% view([0,0])
xlabel('$x_1$ [m]')
ylabel('$x_2$ [m]')
zlabel('$x_3$ [m]')

%%
F_lp = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_lp(:,:,i) = res_lp(end).Q(:,:,i)^(1/2) * S.x;
end


%%
%%% Funnel (SOS Program)
args.xChar = xChar;
args.xN = xNominal;
args.max_iter = 10;
% args.rho = 5; % rho = 0
args.rho = 3; % rho = -1
args.ftol = 5e-4;
args.plot_cost = true;

tic;
[res_sos, cost_sos, rate_sos] = funnel_nonlinear_sos(sys, t, Q0, args);
ctime_all_sos = toc;

F_sos = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_sos(:,:,i) = res_sos(end).step2(:,:,i)^(1/2) * S.x;
end

%%
% %%% HJB equation
% % grid for HJB equation
% nGrid = [151, 151, 151];
% minGrid = [-0.1, -0.1, -0.1];
% maxGrid = [0.1, 0.1, 0.1];
% gr = createGrid(minGrid, maxGrid, nGrid);
% V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) + gr.xs{3}.*gr.xs{3}/Q0(3,3)- 1;
% X0 = Utils.get_level_set(gr, V0, 0.0);
% 
% % solve
% hjb_equation = HJBequation(sys, gr);
% V = zeros([size(V0), length(t)]);
% V(:,:,:,1) = V0;
% tic
% for i = 1:length(t)-1
%     V(:,:,:,i+1) = hjb_equation.solve(V(:,:,:,i), t(i), t(i+1));
% end
% ctime_all_hjb = toc;
% 
% % extract zero-level set
% X = cell(1,length(t));
% X{1} = X0;
% for i = 2:length(t)
%     X{i} = Utils.get_level_set(gr, V(:,:,:,i), 0.0);
% end
%%

%%%%% comparison
figure;
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.5);
h2 = Utils.plot_set(F_lp, t, length(t), 'r', 0.3);
for i = round(linspace(1,length(t),7))
%     tmp = F_sos(:,:,i);
%     h3 = plot3(tmp(1,:), t(i)*ones(1,size(tmp,2)), tmp(2,:), 'color', [76,187,23]/255, 'linewidth', 2);
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
% legend([h1,h2,h3],...
%     '$\mathcal{X}(t)$',....
%     '$\mathcal{F}(t)$ (Proposed)',...
%     '$\mathcal{F}(t)$ (SOS program)',...
%     'location', 'northeast')
title('$\textbf{Comparision of funnels and reachable set}$')

%%
figure;
cla; hold on; grid on;
for i = round(linspace(1,length(t),7))
%     tmp = sys.xN(1:2,i) + F_ltv(1:2,:,i);
%     K = convhull(tmp(1,:)', tmp(2,:)');
%     plot(tmp(1,K), tmp(2,K), 'g', 'linewidth', 2)
    
    tmp = sys.xN(1:2,i) + F_lp(1:2,:,i);
    K = convhull(tmp(1,:)', tmp(2,:)');
    plot(tmp(1,K), tmp(2,K), 'r', 'linewidth', 2)
    
    tmp = sys.xN(1:2,i) + F_sos(1:2,:,i);
    K = convhull(tmp(1,:)', tmp(2,:)');
    plot(tmp(1,K), tmp(2,K), 'b', 'linewidth', 2)
    
    tmp = sys.xN(1:2,i) + xChar(1:2,:,i);
    K = convhull(tmp(1,:)', tmp(2,:)');
    plot(tmp(1,K), tmp(2,K), 'k', 'linewidth', 2)
    
    tmp = sys.xN(1:2,i) + X{i}(1:2,:);
    K = convhull(tmp(1,:)', tmp(2,:)');
    plot(tmp(1,K), tmp(2,K), 'm.', 'linewidth', 2)
    
%     tmp = sys.xN(:,i) + F_lp(:,:,i);
%     plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2)
end
%%
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
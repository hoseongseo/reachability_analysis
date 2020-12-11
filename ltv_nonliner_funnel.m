clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

% given initial set
Q0 = diag([0.25, 0.25]);

% Q0 = diag([0.05; 0.05])^2;

% LTV system
t = linspace(0.0,2.0,51);
A = @(t) [-0.8*t + 0.5, cos(1.5*t + 2); 0.5*t^(2/3), -2*exp(-0.7*t)];
D = @(t) [0.4*cos(t), -0.4*t^2; 0.08*t,  2.8*cos(3*t)];
sys = Dynamics.LTV(2, 0, 2, A, D); 

% Q_ltv = funnel_ltv(sys, linspace(0,0.4,41), Q0);

% wMax = 0.1;
% sys0 = Dynamics.LotkaVolterraNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
% sys = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin


%%%% HJB equation
% grid for HJB equation
nGrid = [101, 101];
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

%%
% args = struct('max_iter', 1, 'ftol', 5e-5, 'plot_cost', true, 'frs', {X}, 'grid', gr);
args = struct('max_iter', 30, 'ftol', 5e-4, 'plot_cost', true, 'grid', gr);
if exist('X')
    args.frs = X;
end
[res, cost, rate] = funnel_nonlinear_lp3(sys, t, Q0, [4,4], [4,4,1,1], [0,0], args);

%%
Q_ltv = funnel_ltv(sys, t, Q0);

S = Utils.Sphere(1,200);
F_ltv = zeros(2, size(S.x,2), length(t));
for i = 1:length(t)
    F_ltv(:,:,i) = Q_ltv(:,:,i)^(1/2) * S.x;
end

S = Utils.Sphere(1,200);
F_lp = zeros(2, size(S.x,2), length(t));
for i = 1:length(t)
    F_lp(:,:,i) = res(end).Q(:,:,i)^(1/2) * S.x;
end


figure(3)
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.8);
h2 = Utils.plot_set(F_ltv, t, length(t), 'r', 0.5);
h3 = Utils.plot_set(F_lp, t, length(t), 'g', 0.3);..., '-', 0.1, 'b');
view([128,11])
camlight left
camlight right
xlabel('$x_1$')
ylabel('$t$ [s]')
zlabel('$x_2$')
xlim([-2.5,2.5])
zlim([-2.5,2.5])
ylim([0,2])
ax = gca;
ax.XLabel.Position = [-0.2194, 2.1422, -2.7370];
ax.YLabel.Position = [2.7314, 0.9289, -2.8754];
legend([h1,h2,h3],...
    '$\underline{\mathcal{X}}(t)$',...
    '$\mathcal{E}(Q_x(t))$',...
    '$\mathcal{F}(t)$',...
    'location', 'northwest')

%%
figure(4)
subplot(2,1,1)
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.8);
h2 = Utils.plot_set(F_ltv, t, length(t), 'r', 0.5);
h3 = Utils.plot_set(F_lp, t, length(t), 'g', 0.3);..., '-', 0.1, 'b');
view([128,11])
camlight left
camlight right
view([90,90])
xlabel('$x_1$')
zlabel('$x_2$')
xlim([-2.5,2.5])
zlim([-2.5,2.5])
ylim([0,2])

subplot(2,1,2)
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.8);
h2 = Utils.plot_set(F_ltv, t, length(t), 'r', 0.5);
h3 = Utils.plot_set(F_lp, t, length(t), 'g', 0.3);..., '-', 0.1, 'b');
view([128,11])
camlight left
camlight right
view([90,0])
xlabel('$x_1$')
ylabel('$t$ [s]')
zlabel('$x_2$')
xlim([-2.5,2.5])
zlim([-2.5,2.5])
ylim([0,2])
% title('$\textbf{Comparision of the funnel and reachable set}$')
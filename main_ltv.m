clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))

set(groot, 'DefaultFigureColor', 'w')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
set(groot, 'DefaultAxesFontSize', 16)

% given initial set
Q0 = diag([0.25, 0.25]);

% LTV system
t = linspace(0,2,401);
A = @(t) [-0.8*t + 0.5, cos(1.5*t + 2); 0.5*t^(2/3), -2*exp(-0.7*t)];
D = @(t) [0.4*cos(t), -0.4*t^2; 0.08*t,  2.8*cos(3*t)];
sys = Dynamics.LTV(2, 0, 2, A, D); 


%%%% Funnel of LTV system (proposed)
Q_ltv = funnel_ltv(sys, t, Q0);

% extract zero level set
S = Utils.Sphere(1,200);
F = zeros(2, size(S.x,2), length(t));
for i = 1:length(t)
    F(:,:,i) = Q_ltv(:,:,i)^(1/2) * S.x;
end

%%%% HJB equation
% grid for HJB equation
nGrid = [201, 201];
minGrid = [-2.5, -2.5];
maxGrid = [2.5, 2.5];
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

%%
%%%%% comparison
figure(3)
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.8);
h2 = Utils.plot_set(F, t, length(t), 'r', 0.5);
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
legend([h1,h2], '$\underline{\mathcal{X}}(t)$', '$\mathcal{E}(Q_x(t))$', 'location', 'northwest')
% title('$\textbf{Comparision of the funnel and reachable set}$')

%%
figure(4)
subplot(2,1,1)
cla; hold on; grid on;
h1 = Utils.plot_set(X, t, length(t), 'k', 0.8);
h2 = Utils.plot_set(F, t, length(t), 'r', 0.5);
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
h2 = Utils.plot_set(F, t, length(t), 'r', 0.5);
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
%%
%%% remove added path
rmpath(genpath('3rd_party/helperOC-master'))
rmpath(genpath('3rd_party/ToolboxLS'))
clear; close all; clc;

addpath(genpath('helperOC-master')) % HJB equation solver
addpath(genpath('ToolboxLS'))


% given initial set
Q0 = diag([0.25, 0.25]);

% LTV system
t = linspace(0,2,201);
A = @(t) [-0.8*t + 0.5, cos(1.5*t + 2); 0.5*t^(2/3), -2*exp(-0.7*t)];
D = @(t) [0.4*cos(t), -0.4*t^2; 0.08*t,  2.8*cos(3*t)];
sys = Dynamics.LTV(2, 0, 2, A, D); 


%%%% Funnel of LTV system (proposed)
funnel = FunnelLTV(sys, Q0, t);

% extract zero level set
S = Utils.Sphere(1,200);
F = zeros(2, size(S.x,2), length(t));
for i = 1:length(t)
    F(:,:,i) = funnel.Q(:,:,i)^(1/2) * S.x;
end

%%%% HJB equation
% grid for HJB equation
nGrid = [101, 101];
minGrid = [-2.5, -2.5];
maxGrid = [2.5, 2.5];
gr = createGrid(minGrid, maxGrid, nGrid);
V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) - 1;
X0 = getLevelSet(gr, V0, 0.0);

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
    X{i} = getLevelSet(gr, V(:,:,i), 0.0);
end

%%%%% comparison
figure;
cla; hold on; grid on;
h1 = plotSet(X, t, length(t), 'k', 0.7);
h2 = plotSet(F, t, length(t), 'r', 0.5);
view([128,11])
camlight left
camlight right
xlabel('$x_1$')
ylabel('$t$ [s]')
zlabel('$x_2$')
ax = gca;
ax.XLabel.Position = [-0.1725, 2.1516, -3.3917];
ax.YLabel.Position = [2.2098, 0.9221, -3.5364];
legend([h1,h2], '$\underline{\mathcal{X}}(t)$', '$\mathcal{E}(Q_x(t))$', 'location', 'northwest')

%%% remove added path
rmpath(genpath('helperOC-master'))
rmpath(genpath('ToolboxLS'))
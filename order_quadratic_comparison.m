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
S = Utils.Sphere(1,1000);

%%% system (polynomial dynamics)
t = linspace(0,1,151);
sys0 = Dynamics.LotkaVolterraNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
sys = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin

%%
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

A = zeros(1,length(t));
for i = 1:length(t)
    A(i) = area( polyshape(X{i}') );
end

%%% Funnel (proposed)
%%
args = struct('max_iter', 3, 'ftol', 5e-5, 'plot_cost', false);
cvx_precision best
% cvx_solver 

tic
[res_quad, cost_quad, rate_quad] = funnel_nonlinear_quad(sys, t, Q0, args);
ctime_quad = toc;
F_quad = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_quad(:,:,i) = res_quad(end).Q(:,:,i)^(1/2) * S.x;
end
A_quad = zeros(1,length(t));
for i = 1:length(t)
    A_quad(i) = area( polyshape(F_quad(:,:,i)') );
end

%%
tic
[res_22, cost_22, rate_22] = funnel_nonlinear_lp(sys, t, Q0, [2,2], args);
ctime_22 = toc;
F_22 = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_22(:,:,i) = res_22(end).Q(:,:,i)^(1/2) * S.x;
end
A_22 = zeros(1,length(t));
for i = 1:length(t)
    A_22(i) = area( polyshape(F_22(:,:,i)') );
end


%%
tic;
[res_33, cost_33, rate_33] = funnel_nonlinear_lp(sys, t, Q0, [3,3], args);
ctime_33 = toc;
F_33 = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_33(:,:,i) = res_33(end).Q(:,:,i)^(1/2) * S.x;
end
A_33 = zeros(1,length(t));
for i = 1:length(t)
    A_33(i) = area( polyshape(F_33(:,:,i)') );
end


%%
tic;
[res_44, cost_44, rate_44] = funnel_nonlinear_lp(sys, t, Q0, [4,4], args);
ctime_44 = toc;
F_44 = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_44(:,:,i) = res_44(end).Q(:,:,i)^(1/2) * S.x;
end
A_44 = zeros(1,length(t));
for i = 1:length(t)
    A_44(i) = area( polyshape(F_44(:,:,i)') );
end

%%
tic;
[res_55, cost_55, rate_55] = funnel_nonlinear_lp(sys, t, Q0, [5,5], args);
ctime_55 = toc;
F_55 = zeros([size(S.x), length(t)]);
for i = 1:length(t)
    F_55(:,:,i) = res_55(end).Q(:,:,i)^(1/2) * S.x;
end
A_55 = zeros(1,length(t));
for i = 1:length(t)
    A_55(i) = area( polyshape(F_55(:,:,i)') );
end

%%

figure(1)
cla; hold on; grid on;
% plot(t, A, 'k', 'linewidth', 2)
h1 = plot(t, A_22./A, 'b', 'linewidth', 2);
h2 = plot(t, A_33./A, 'r', 'linewidth', 2);
h3 = plot(t, A_44./A, 'g', 'linewidth', 2);
h4 = plot(t, A_55./A, 'm', 'linewidth', 2);
axis tight;
xlabel('$t$ [s]')
ylabel('Conservativeness')
% title('Vol$(\mathcal{F}(t))/$Vol$(\mathcal{X}(t))$')
legend([h1,h2,h3,h4],...
    '$N = [2,2]$',...
    '$N = [3,3]$',...
    '$N = [4,4]$',...
    '$N = [5,5]$',...
    'location', 'northwest')
%%
figure(3)
cla; hold on; grid on;
yyaxis left
plot([ctime_22, ctime_33, ctime_44, ctime_55], 'b-o', 'linewidth', 2)
yyaxis right
plot([sum(A_22./A)*(t(2)-t(1)), sum(A_33./A)*(t(2)-t(1)), sum(A_44./A)*(t(2)-t(1)), sum(A_55./A)*(t(2)-t(1))], 'r-o', 'linewidth', 2)
ax = gca;
ax.XTick = [1.0, 2.0, 3.0, 4.0];
ax.XTickLabel = {'$N=[2,2]$', '$N=[3,3]$', '$N=[4,4]$', '$N=[5,5]$'};
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
ax.YAxis(1).Label.String = 'Computation time [s]';
ax.YAxis(2).Label.String = 'Conservativeness';
ax.Position = [0.11, 0.11, 0.7665, 0.8150];
%%
bdry_X = sys.xN(:,1) + X{1};
for i = 2:length(t)
    tmp = [bdry_X, sys.xN(:,i) + X{i}];
    K = boundary(tmp(1,:)', tmp(2,:)', 0.999);
    bdry_X = tmp(:,K);
    disp(['i = ',num2str(i)])
end

%%
bdry_22 = sys.xN(:,1) + F_22(:,:,1);
for i = 2:length(t)
    tmp = [bdry_22, sys.xN(:,i) + F_22(:,:,i)];
    K = boundary(tmp(1,:)', tmp(2,:)', 0.999);
    bdry_22 = tmp(:,K);
    disp(['i = ',num2str(i)])
end
%%
bdry_33 = sys.xN(:,1) + F_33(:,:,1);
for i = 2:length(t)
    tmp = [bdry_33, sys.xN(:,i) + F_33(:,:,i)];
    K = boundary(tmp(1,:)', tmp(2,:)', 0.995);
    bdry_33 = tmp(:,K);
    disp(['i = ',num2str(i)])
end
%%
bdry_44 = sys.xN(:,1) + F_44(:,:,1);
for i = 2:length(t)
    tmp = [bdry_44, sys.xN(:,i) + F_44(:,:,i)];
    K = boundary(tmp(1,:)', tmp(2,:)', 0.995);
    bdry_44 = tmp(:,K);
    disp(['i = ',num2str(i)])
end
%%
bdry_55 = sys.xN(:,1) + F_55(:,:,1);
for i = 2:length(t)
    tmp = [bdry_55, sys.xN(:,i) + F_55(:,:,i)];
    K = boundary(tmp(1,:)', tmp(2,:)', 0.993);
    bdry_55 = tmp(:,K);
    disp(['i = ',num2str(i)])
end
%%
figure(2)
cla; hold on; grid on; axis equal
% plot(bdry_X(1,:), bdry_X(2,:), 'k', 'linewidth', 2)
h1 = patch('xdata', bdry_22(1,:), 'ydata', bdry_22(2,:), 'linestyle', '-', 'edgecolor', 'r', 'linewidth', 0.1, 'facecolor', 'r', 'facealpha', 0.5);
% patch('xdata', bdry_33(1,:), 'ydata', bdry_33(2,:), 'linestyle', '-', 'edgecolor', 'r', 'linewidth', 1, 'facecolor', 'r', 'facealpha', 0.2)
% patch('xdata', bdry_44(1,:), 'ydata', bdry_44(2,:), 'linestyle', '-', 'edgecolor', [255,129,0]/255, 'linewidth', 1, 'facecolor', [255,129,0]/255, 'facealpha', 0.2)
h2 = patch('xdata', bdry_55(1,:), 'ydata', bdry_55(2,:), 'linestyle', '-', 'edgecolor', [76,187,23]/255, 'linewidth', 0.1, 'facecolor', [76,187,23]/255, 'facealpha', 0.5);
h3 = patch('xdata', bdry_X(1,:), 'ydata', bdry_X(2,:),   'linestyle', '-', 'edgecolor', 'k', 'linewidth', 0.1, 'facecolor', 'k', 'facealpha', 0.5);

for i = round(linspace(1,length(t),5))
    tmp = sys.xN(:,i) + F_22(:,:,i);
    plot(tmp(1,:), tmp(2,:), 'r', 'linewidth', 2);
    
%     tmp = sys.xN(:,i) + F_33(:,:,i);
%     h2 = plot(tmp(1,:), tmp(2,:), 'r', 'linewidth', 2);
    
%     tmp = sys.xN(:,i) + F_44(:,:,i);
%     h3 = plot(tmp(1,:), tmp(2,:), 'color', [255,129,0]/255, 'linewidth', 2);
    
    tmp = sys.xN(:,i) + X{i};
    plot(tmp(1,:), tmp(2,:), 'k', 'linewidth', 1);
    
%     tmp = sys.xN(:,i) + F_55(:,:,i);
%     plot(tmp(1,:), tmp(2,:), 'color', [76,187,23]/255, 'linestyle', '-', 'linewidth', 2);
    
    
    text(sys.xN(1,i), sys.xN(2,i), ['$t = ',num2str(round(t(i),2)),'$'],...
        'horizontalalignment', 'center', 'fontsize', 16, 'color', 'w')
end

xlabel('$x_1$')
ylabel('$x_2$')
legend([h3,h1,h2],...
    '$\mathcal{X}(t)$',...
    '$\mathcal{F}(t)$, $N=[2,2]$',...
    '$\mathcal{F}(t)$, $N=[5,5]$',...
    'location', 'southeast')
    
%% Comparison with quadratic approximation
figure(2)
cla; hold on; grid on; axis equal
% plot(bdry_X(1,:), bdry_X(2,:), 'k', 'linewidth', 2)
% h1 = patch('xdata', bdry_22(1,:), 'ydata', bdry_22(2,:), 'linestyle', '-', 'edgecolor', 'r', 'linewidth', 0.1, 'facecolor', 'r', 'facealpha', 0.5);
% patch('xdata', bdry_33(1,:), 'ydata', bdry_33(2,:), 'linestyle', '-', 'edgecolor', 'r', 'linewidth', 1, 'facecolor', 'r', 'facealpha', 0.2)
% patch('xdata', bdry_44(1,:), 'ydata', bdry_44(2,:), 'linestyle', '-', 'edgecolor', [255,129,0]/255, 'linewidth', 1, 'facecolor', [255,129,0]/255, 'facealpha', 0.2)
% h2 = patch('xdata', bdry_55(1,:), 'ydata', bdry_55(2,:), 'linestyle', '-', 'edgecolor', [76,187,23]/255, 'linewidth', 0.1, 'facecolor', [76,187,23]/255, 'facealpha', 0.5);
% h3 = patch('xdata', bdry_X(1,:), 'ydata', bdry_X(2,:),   'linestyle', '-', 'edgecolor', 'k', 'linewidth', 0.1, 'facecolor', 'k', 'facealpha', 0.5);

for i = round(linspace(1,length(t),6))
    tmp = sys.xN(:,i) + F_quad(:,:,i);
    h1 = plot(tmp(1,:), tmp(2,:), 'b', 'linewidth', 2);
    
    tmp = sys.xN(:,i) + F_22(:,:,i);
    h2 = plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2);
    
%     tmp = sys.xN(:,i) + F_33(:,:,i);
%     h2 = plot(tmp(1,:), tmp(2,:), 'r', 'linewidth', 2);
    
%     tmp = sys.xN(:,i) + F_44(:,:,i);
%     h3 = plot(tmp(1,:), tmp(2,:), 'color', [255,129,0]/255, 'linewidth', 2);
    
    tmp = sys.xN(:,i) + X{i};
    h3 = plot(tmp(1,:), tmp(2,:), 'k', 'linewidth', 2);
    
    text(sys.xN(1,i), sys.xN(2,i), ['$t = ',num2str(round(t(i),2)),'$'],...
        'horizontalalignment', 'center', 'fontsize', 16, 'color', 'k')
end

xlabel('$x_1$')
ylabel('$x_2$')
legend([h3,h1,h2],...
    '$\mathcal{X}(t)$',...
    '$\mathcal{F}(t)$, (quadratic function)',...
    '$\mathcal{F}(t)$, $N=[2,2]$',...
    'location', 'southeast')
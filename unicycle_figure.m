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
Q0 = (diag([0.03,0.03,0.03])^2);

% disturbance bound
wMax = 0.05;

% sphere
S = Utils.Sphere(2,500);

%%% system (polynomial dynamics)
t = linspace(0,1,51);
% t = linspace(0,0.1,11);


%% case 1
rho = -1;
u = 3*rho*sin(2*pi*t);
sys0 = Dynamics.UnicycleNominal([0; 0; 0], u, t); % nominal dynamics with the initial condition
sys = Dynamics.Unicycle(sys0, wMax); % system shifted to the origin

load('unicycle_test2_rm3.mat')
load('res_sos_rm3.mat')
xNominal_m3 = sys.xN + [0; -0.14; 0];
xChar_m3 = xChar;
tic
res_lp = funnel_nonlinear_lp4(sys, t, Q0, [2,2,2], [0,0,0,0], [2,2,2]);
toc
Q_lp_m3 = res_lp(end).Q;
Q_sos_m3 = res_sos(end).step2;
Q_sos_3 = res_sos(end).step2;

%% case 2
rho = -0.5;
u = 3*rho*sin(2*pi*t);
sys0 = Dynamics.UnicycleNominal([0; 0; 0], u, t); % nominal dynamics with the initial condition
sys = Dynamics.Unicycle(sys0, wMax); % system shifted to the origin

load('unicycle_test2_rm15.mat')
load('res_sos_rm15.mat')
xNominal_m15 = sys.xN + [0; -0.07; 0];
xChar_m15 = xChar;

res_lp = funnel_nonlinear_lp4(sys, t, Q0, [2,2,2], [0,0,0,0], [2,2,2]);
Q_lp_m15 = res_lp(end).Q;
Q_sos_m15 = res_sos(end).step2;
Q_sos_15 = res_sos(end).step2;

%% case 3
rho = 0;
u = 3*rho*sin(2*pi*t);
sys0 = Dynamics.UnicycleNominal([0; 0; 0], u, t); % nominal dynamics with the initial condition
sys = Dynamics.Unicycle(sys0, wMax); % system shifted to the origin

load('unicycle_test2_r0.mat')
load('res_sos_r0.mat')
xNominal_0 = sys.xN;
xChar_0 = xChar;

res_lp = funnel_nonlinear_lp4(sys, t, Q0, [2,2,2], [0,0,0,0], [2,2,2]);
Q_lp_0 = res_lp(end).Q;
Q_sos_0 = res_sos(end).step2;

%% case 4
rho = 0.5;
u = 3*rho*sin(2*pi*t);
sys0 = Dynamics.UnicycleNominal([0; 0; 0], u, t); % nominal dynamics with the initial condition
sys = Dynamics.Unicycle(sys0, wMax); % system shifted to the origin

load('unicycle_test2_r15.mat')
xNominal_15 = sys.xN + [0; 0.07; 0];
xChar_15 = xChar;

res_lp = funnel_nonlinear_lp4(sys, t, Q0, [2,2,2], [0,0,0,0], [2,2,2]);
Q_lp_15 = res_lp(end).Q;

%% case 5
rho = 1.0;
u = 3*rho*sin(2*pi*t);
sys0 = Dynamics.UnicycleNominal([0; 0; 0], u, t); % nominal dynamics with the initial condition
sys = Dynamics.Unicycle(sys0, wMax); % system shifted to the origin

load('unicycle_test2_r3.mat')
xNominal_3 = sys.xN  + [0; 0.14; 0];
xChar_3 = xChar;

res_lp = funnel_nonlinear_lp4(sys, t, Q0, [2,2,2], [0,0,0,0], [2,2,2]);
Q_lp_3 = res_lp(end).Q;


%%
figure(1)
cla; hold on; grid on; 
axis equal; 
axis tight;
% tmp0 = xNominal_0;
% tmp = xChar_0;
% for k = round(linspace(1,size(xChar,2),101))%1:size(xChar,2)
%     plot(tmp0(2,:) + squeeze(tmp(2,k,:))', tmp0(1,:) + squeeze(tmp(1,k,:))', 'color', 0.3*ones(3,1));
% end
% 
% tmp0 = xNominal_15;
% tmp = xChar_15;
% for k = round(linspace(1,size(xChar,2),101))%1:size(xChar,2)
%     plot(tmp0(2,:) + squeeze(tmp(2,k,:))', tmp0(1,:) + squeeze(tmp(1,k,:))', 'color', 0.3*ones(3,1));
% end
% 
% tmp0 = xNominal_3;
% tmp = xChar_3;
% for k = round(linspace(1,size(xChar,2),101))%1:size(xChar,2)
%     plot(tmp0(2,:) + squeeze(tmp(2,k,:))', tmp0(1,:) + squeeze(tmp(1,k,:))', 'color', 0.3*ones(3,1));
% end
% 
% tmp0 = xNominal_m15;
% tmp = xChar_m15;
% for k = round(linspace(1,size(xChar,2),101))%1:size(xChar,2)
%     plot(tmp0(2,:) + squeeze(tmp(2,k,:))', tmp0(1,:) + squeeze(tmp(1,k,:))', 'color', 0.3*ones(3,1));
% end
% 
% tmp0 = xNominal_m3;
% tmp = xChar_m3;
% for k = round(linspace(1,size(xChar,2),101))%1:size(xChar,2)
%     plot(tmp0(2,:) + squeeze(tmp(2,k,:))', tmp0(1,:) + squeeze(tmp(1,k,:))', 'color', 0.3*ones(3,1));
% end


tmp0 = xNominal_0;
tmpQ_lp = Q_lp_0;
tmpQ_sos = Q_sos_0;
tmp = xChar_0;
X0_lp_all = [];
X0_sos_all = [];
for idx = 1:size(xNominal,2)
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_lp_all = [X0_lp_all, X0_lp];
%     X0_sos_all = [X0_sos_all, X0_sos];
    
    if idx > 1
        K = boundary(X0_lp_all(1,:)', X0_lp_all(2,:)', 0.95);
        X0_lp_all = X0_lp_all(:,K);
%         K = boundary(X0_sos_all(1,:)', X0_sos_all(2,:)', 0.98);
%         X0_sos_all = X0_sos_all(:,K);
    end
    disp(['k = ', num2str(idx)])
end
h1 = patch('xdata', X0_lp_all(2,:), 'ydata', X0_lp_all(1,:),...
    'facecolor', 'r', 'facealpha', 0.2,...
    'linestyle', 'none', 'edgecolor', 'r', 'linewidth', 2);
% patch('xdata', X0_sos_all(2,:), 'ydata', X0_sos_all(1,:),...
%     'facecolor', [76,187,23]/255, 'facealpha', 0.0,...
%     'linestyle', '-', 'edgecolor', [76,187,23]/255, 'linewidth', 2)
for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    patch('xdata', X0_lp(2,:), 'ydata', X0_lp(1,:),...
        'facecolor', 'r', 'facealpha', 0.0,...
        'linestyle', '-', 'edgecolor', 'r', 'linewidth', 2)
    X0_sos = tmp0(1:2,idx) + tmpQ_sos(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
%     h2 = patch('xdata', X0_sos(2,:), 'ydata', X0_sos(1,:),...
%         'facecolor', [0,151,50]/255, 'facealpha', 0.0,...
%         'linestyle', '-', 'edgecolor', [0,151,50]/255, 'linewidth', 2);
    h2 = plot(X0_sos(2,:), X0_sos(1,:), 'color', [0,151,50]/255, 'linewidth', 2);
end

for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
    K = convhull(tmp(1,:,idx)', tmp(2,:,idx)');
    h3 = patch('xdata', tmp0(2,idx) + tmp(2,K,idx), 'ydata', tmp0(1,idx) + tmp(1,K,idx),...
        'facecolor', 'k', 'facealpha', 0.5,...
        'linestyle', 'none', 'edgecolor', 'k', 'linewidth', 1);
%     plot(tmp0(2,idx) + tmp(2,K,idx),...
%          tmp0(1,idx) + tmp(1,K,idx), '-', 'color', 'k')
    plot(tmp0(2,idx) + tmp(2,K,idx), tmp0(1,idx) + tmp(1,K,idx), 'k', 'linewidth', 1)
end


tmp0 = xNominal_15;
tmpQ_lp = Q_lp_15;
tmpQ_sos = Q_sos_15;
tmp = xChar_15;
X0_lp_all = [];
X0_sos_all = [];
for idx = 1:size(xNominal,2)
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_lp_all = [X0_lp_all, X0_lp];
%     X0_sos_all = [X0_sos_all, X0_sos];
    
    if idx > 1
        K = boundary(X0_lp_all(1,:)', X0_lp_all(2,:)', 0.95);
        X0_lp_all = X0_lp_all(:,K);
%         K = boundary(X0_sos_all(1,:)', X0_sos_all(2,:)', 0.98);
%         X0_sos_all = X0_sos_all(:,K);
    end
    disp(['k = ', num2str(idx)])
end
patch('xdata', X0_lp_all(2,:), 'ydata', X0_lp_all(1,:),...
    'facecolor', 'r', 'facealpha', 0.2,...
    'linestyle', 'none', 'edgecolor', 'r', 'linewidth', 2)
% patch('xdata', X0_sos_all(2,:), 'ydata', X0_sos_all(1,:),...
%     'facecolor', [76,187,23]/255, 'facealpha', 0.0,...
%     'linestyle', '-', 'edgecolor', [76,187,23]/255, 'linewidth', 2)
for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
        
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    patch('xdata', X0_lp(2,:), 'ydata', X0_lp(1,:),...
        'facecolor', 'r', 'facealpha', 0.0,...
        'linestyle', '-', 'edgecolor', 'r', 'linewidth', 2)
    
%     X0_sos = tmp0(1:2,idx) + tmpQ_sos(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_sos = tmpQ_sos(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_sos = tmp0(1:2,idx) + [X0_sos(1,:); -X0_sos(2,:)];
%     h2 = patch('xdata', X0_sos(2,:), 'ydata', X0_sos(1,:),...
%         'facecolor', [0,151,50]/255, 'facealpha', 0.0,...
%         'linestyle', '-', 'edgecolor', [0,151,50]/255, 'linewidth', 2);
    plot(X0_sos(2,:), X0_sos(1,:), 'color', [0,151,50]/255, 'linewidth', 2);
end

for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
    K = convhull(tmp(1,:,idx)', tmp(2,:,idx)');
    patch('xdata', tmp0(2,idx) + tmp(2,K,idx), 'ydata', tmp0(1,idx) + tmp(1,K,idx),...
        'facecolor', 'k', 'facealpha', 0.5,...
        'linestyle', 'none', 'edgecolor', 'k', 'linewidth', 1)
%     plot(tmp0(2,idx) + tmp(2,K,idx),...
%          tmp0(1,idx) + tmp(1,K,idx), '-', 'color', 'k')
    plot(tmp0(2,idx) + tmp(2,K,idx), tmp0(1,idx) + tmp(1,K,idx), 'k', 'linewidth', 1)
end


tmp0 = xNominal_3;
tmpQ_lp = Q_lp_3;
tmpQ_sos = Q_sos_3;
tmp = xChar_3;
X0_lp_all = [];
X0_sos_all = [];
for idx = 1:size(xNominal,2)
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_lp_all = [X0_lp_all, X0_lp];
%     X0_sos_all = [X0_sos_all, X0_sos];
    
    if idx > 1
        K = boundary(X0_lp_all(1,:)', X0_lp_all(2,:)', 0.95);
        X0_lp_all = X0_lp_all(:,K);
%         K = boundary(X0_sos_all(1,:)', X0_sos_all(2,:)', 0.98);
%         X0_sos_all = X0_sos_all(:,K);
    end
    disp(['k = ', num2str(idx)])
end
patch('xdata', X0_lp_all(2,:), 'ydata', X0_lp_all(1,:),...
    'facecolor', 'r', 'facealpha', 0.2,...
    'linestyle', 'none', 'edgecolor', 'r', 'linewidth', 2)
% patch('xdata', X0_sos_all(2,:), 'ydata', X0_sos_all(1,:),...
%     'facecolor', [76,187,23]/255, 'facealpha', 0.0,...
%     'linestyle', '-', 'edgecolor', [76,187,23]/255, 'linewidth', 2)
for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
        
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    patch('xdata', X0_lp(2,:), 'ydata', X0_lp(1,:),...
        'facecolor', 'r', 'facealpha', 0.0,...
        'linestyle', '-', 'edgecolor', 'r', 'linewidth', 2)
    
%     X0_sos = tmp0(1:2,idx) + tmpQ_sos(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_sos = tmpQ_sos(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_sos = tmp0(1:2,idx) + [X0_sos(1,:); -X0_sos(2,:)];
%     h2 = patch('xdata', X0_sos(2,:), 'ydata', X0_sos(1,:),...
%         'facecolor', [0,151,50]/255, 'facealpha', 0.0,...
%         'linestyle', '-', 'edgecolor', [0,151,50]/255, 'linewidth', 2);
    plot(X0_sos(2,:), X0_sos(1,:), 'color', [0,151,50]/255, 'linewidth', 2);
end

for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
    K = convhull(tmp(1,:,idx)', tmp(2,:,idx)');
    patch('xdata', tmp0(2,idx) + tmp(2,K,idx), 'ydata', tmp0(1,idx) + tmp(1,K,idx),...
        'facecolor', 'k', 'facealpha', 0.5,...
        'linestyle', 'none', 'edgecolor', 'k', 'linewidth', 1)
%     plot(tmp0(2,idx) + tmp(2,K,idx),...
%          tmp0(1,idx) + tmp(1,K,idx), '-', 'color', 'k')
    plot(tmp0(2,idx) + tmp(2,K,idx), tmp0(1,idx) + tmp(1,K,idx), 'k', 'linewidth', 1)
end

tmp0 = xNominal_m15;
tmpQ_lp = Q_lp_m15;
tmpQ_sos = Q_sos_m15;
tmp = xChar_m15;
X0_lp_all = [];
X0_sos_all = [];
for idx = 1:size(xNominal,2)
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_lp_all = [X0_lp_all, X0_lp];
%     X0_sos_all = [X0_sos_all, X0_sos];
    
    if idx > 1
        K = boundary(X0_lp_all(1,:)', X0_lp_all(2,:)', 0.95);
        X0_lp_all = X0_lp_all(:,K);
%         K = boundary(X0_sos_all(1,:)', X0_sos_all(2,:)', 0.98);
%         X0_sos_all = X0_sos_all(:,K);
    end
    disp(['k = ', num2str(idx)])
end
patch('xdata', X0_lp_all(2,:), 'ydata', X0_lp_all(1,:),...
    'facecolor', 'r', 'facealpha', 0.2,...
    'linestyle', 'none', 'edgecolor', 'r', 'linewidth', 2)
% patch('xdata', X0_sos_all(2,:), 'ydata', X0_sos_all(1,:),...
%     'facecolor', [76,187,23]/255, 'facealpha', 0.0,...
%     'linestyle', '-', 'edgecolor', [76,187,23]/255, 'linewidth', 2)
for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
        
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    patch('xdata', X0_lp(2,:), 'ydata', X0_lp(1,:),...
        'facecolor', 'r', 'facealpha', 0.0,...
        'linestyle', '-', 'edgecolor', 'r', 'linewidth', 2)
    
    X0_sos = tmp0(1:2,idx) + tmpQ_sos(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
%     h2 = patch('xdata', X0_sos(2,:), 'ydata', X0_sos(1,:),...
%         'facecolor', [0,151,50]/255, 'facealpha', 0.0,...
%         'linestyle', '-', 'edgecolor', [0,151,50]/255, 'linewidth', 2);
    plot(X0_sos(2,:), X0_sos(1,:), 'color', [0,151,50]/255, 'linewidth', 2);
end
for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
    K = convhull(tmp(1,:,idx)', tmp(2,:,idx)');
    patch('xdata', tmp0(2,idx) + tmp(2,K,idx), 'ydata', tmp0(1,idx) + tmp(1,K,idx),...
        'facecolor', 'k', 'facealpha', 0.5,...
        'linestyle', 'none', 'edgecolor', 'k', 'linewidth', 1)
%     plot(tmp0(2,idx) + tmp(2,K,idx),...
%          tmp0(1,idx) + tmp(1,K,idx), '-', 'color', 'k')
    plot(tmp0(2,idx) + tmp(2,K,idx), tmp0(1,idx) + tmp(1,K,idx), 'k', 'linewidth', 1)
end

tmp0 = xNominal_m3;
tmpQ_lp = Q_lp_m3;
tmpQ_sos = Q_sos_m3;
tmp = xChar_m3;
X0_lp_all = [];
X0_sos_all = [];
for idx = 1:size(xNominal,2)
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    X0_lp_all = [X0_lp_all, X0_lp];
%     X0_sos_all = [X0_sos_all, X0_sos];
    
    if idx > 1
        K = boundary(X0_lp_all(1,:)', X0_lp_all(2,:)', 0.95);
        X0_lp_all = X0_lp_all(:,K);
%         K = boundary(X0_sos_all(1,:)', X0_sos_all(2,:)', 0.98);
%         X0_sos_all = X0_sos_all(:,K);
    end
    disp(['k = ', num2str(idx)])
end
patch('xdata', X0_lp_all(2,:), 'ydata', X0_lp_all(1,:),...
    'facecolor', 'r', 'facealpha', 0.2,...
    'linestyle', 'none', 'edgecolor', 'r', 'linewidth', 2)
% patch('xdata', X0_sos_all(2,:), 'ydata', X0_sos_all(1,:),...
%     'facecolor', [76,187,23]/255, 'facealpha', 0.0,...
%     'linestyle', '-', 'edgecolor', [76,187,23]/255, 'linewidth', 2)
for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
        
    X0_lp = tmp0(1:2,idx) + tmpQ_lp(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
    patch('xdata', X0_lp(2,:), 'ydata', X0_lp(1,:),...
        'facecolor', 'r', 'facealpha', 0.0,...
        'linestyle', '-', 'edgecolor', 'r', 'linewidth', 2)
    
    X0_sos = tmp0(1:2,idx) + tmpQ_sos(1:2,1:2,idx)^(1/2) * Math.Sphere(1,200).x;
%     h2 = patch('xdata', X0_sos(2,:), 'ydata', X0_sos(1,:),...
%         'facecolor', [0,151,50]/255, 'facealpha', 0.0,...
%         'linestyle', '-', 'edgecolor', [0,151,50]/255, 'linewidth', 2);
    plot(X0_sos(2,:), X0_sos(1,:), 'color', [0,151,50]/255, 'linewidth', 2);
end

for idx = round(linspace(1,size(xNominal,2),11))%1:size(xNominal,2)%round(linspace(1,size(xNominal,2),21))
    K = convhull(tmp(1,:,idx)', tmp(2,:,idx)');
    patch('xdata', tmp0(2,idx) + tmp(2,K,idx), 'ydata', tmp0(1,idx) + tmp(1,K,idx),...
        'facecolor', 'k', 'facealpha', 0.5,...
        'linestyle', 'none', 'edgecolor', 'k', 'linewidth', 1)
%     plot(tmp0(2,idx) + tmp(2,K,idx),...
%          tmp0(1,idx) + tmp(1,K,idx), '-', 'color', 'k')
    plot(tmp0(2,idx) + tmp(2,K,idx), tmp0(1,idx) + tmp(1,K,idx), 'k', 'linewidth', 1)
end
lgd = legend([h3,h2,h1],...
    '$\mathcal{X}(t)$ (HJB equation)',...
    '$\mathcal{F}(t)$ (SOS program)',...
    '$\mathcal{F}(t)$ (Proposed)',...
    'location', 'southeast');
lgd.Position = [0.6447, 0.1298, 0.3089, 0.1353];

xlabel('$x_2$ [m]')
ylabel('$x_1$ [m]')


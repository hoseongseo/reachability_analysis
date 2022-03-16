clear; close all; clc;

t = linspace(0,1,101);
sys0 = Dynamics.LotkaVolterraContNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
cmap = winter;


nPts = 20;
nPts2 = 200;

% cmap = [linspace(1,0.5,5)', linspace(0,0.5,5)', linspace(0,0.5,5)'];
figure(11)
cla; hold on; grid on; axis equal; axis tight;
set(gcf, 'position', [562         569        1353         532])

load('FRS.mat')
for k = 1:20:length(X)
    %%% HJB eqn
    x_ = sys0.xN(:,k) + X{k};
    h0 = plot(x_(1,:), x_(2,:), 'k', 'linewidth', 2);
    
%     text(sys0.xN(1,k), sys0.xN(2,k), ['$t=',num2str(t(k)),'$ s'],...
%         'horizontalalignment', 'center', 'fontsize', 14)
end

load('discrete_100.mat')
idx = 1:20:size(Q_proposed,3);
% col = cmap(5,:);
col = 'c';
% col = Utils.getColor(5, [1,5], cmap);
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h5 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

load('discrete_50.mat')
idx = 1:10:size(Q_proposed,3);
% col = cmap(4,:);
col = 'm';
% col = Utils.getColor(4, [1,5], cmap);
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*2+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h4 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*2+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

load('discrete_20.mat')
idx = 1:4:size(Q_proposed,3);
% col = cmap(3,:);
col = 'b';
% col = Utils.getColor(3, [1,5], cmap);
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*5+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h3 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*5+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

load('discrete_10.mat')
idx = 1:2:size(Q_proposed,3);
% col = cmap(2,:);
col = 'g';
% col = Utils.getColor(2, [1,5], cmap);
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*10+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h2 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*10+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

load('discrete_5.mat')
idx = 1:size(Q_proposed,3);
% col = cmap(1,:);
col = 'r';
% col = Utils.getColor(1, [1,5], cmap);
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*20+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h1 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*20+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

k = 1;
text(sys0.xN(1,k), sys0.xN(2,k), ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 21;
text(sys0.xN(1,k), sys0.xN(2,k)-0.09, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 41;
text(sys0.xN(1,k), sys0.xN(2,k)-0.1, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 61;
text(sys0.xN(1,k), sys0.xN(2,k)-0.12, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 81;
text(sys0.xN(1,k), sys0.xN(2,k)-0.13, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 101;
text(sys0.xN(1,k), sys0.xN(2,k)-0.15, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
    
legend([h0,h1,h2,h3,h4,h5],...
    'HJB PDE',...
    '$N = 5$',...
    '$N = 10$',...
    '$N = 20$',...
    '$N = 50$',...
    '$N = 100$',...
    'location', 'northeast')

xlabel('$x_1$')
ylabel('$x_2$')

set(gca, 'xlim', [0.05, 1.25])
set(gca, 'ylim', [1.05, 1.6])

%%
figure(12)
cla; hold on; grid on; axis equal; axis tight;
set(gcf, 'position', [562         569        1353         532])

t = linspace(0,1,101);
sys0 = Dynamics.LotkaVolterraContNominal([1.2; 1.1], t); % nominal dynamics with the initial condition


load('FRS.mat')
for k = 1:20:length(X)
    %%% HJB eqn
    x_ = sys0.xN(:,k) + X{k};
    h0 = plot(x_(1,:), x_(2,:), 'k', 'linewidth', 2);
    
%     text(q_proposed(1,k), q_proposed(2,k), ['$t=',num2str(t(k)),'$ s'],...
%         'horizontalalignment', 'center', 'fontsize', 14)
end

load('cont_5.mat')
idx = 1:size(Q_proposed,3);
% col = cmap(1,:);
col = 'r';
% col = Utils.getColor(1, [1,5], cmap);
% t = linspace(0,1,size(Q_proposed,3));
% sys0 = Dynamics.LotkaVolterraContNominal([1.2; 1.1], t); % nominal dynamics with the initial condition

for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*20+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h1 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*20+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

load('cont_10.mat')
idx = 1:2:size(Q_proposed,3);
% col = cmap(2,:);
col = 'g';
% col = Utils.getColor(2, [1,5], cmap);
% t = linspace(0,1,size(Q_proposed,3));
% sys0 = Dynamics.LotkaVolterraContNominal([1.2; 1.1], t); % nominal dynamics with the initial condition

for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*10+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h2 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*10+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

load('cont_20.mat')
idx = 1:4:size(Q_proposed,3);
% col = cmap(3,:)
col = 'b';
% col = Utils.getColor(3, [1,5], cmap);
% t = linspace(0,1,size(Q_proposed,3));
% sys0 = Dynamics.LotkaVolterraContNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*5+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h3 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*5+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
    
%     %%% HJB eqn
%     x_ = sys.xN(:,k) + X{k};
%     h6 = plot(x_(1,:), x_(2,:), 'k-', 'linewidth', 2);
    
%     text(q_proposed(1,k), q_proposed(2,k), ['$t=',num2str(t(k)),'$ s'],...
%         'horizontalalignment', 'center', 'fontsize', 14)
end

load('cont_50.mat')
idx = 1:10:size(Q_proposed,3);
% col = cmap(4,:);
col = 'm';
% col = Utils.getColor(4, [1,5], cmap);
% t = linspace(0,1,size(Q_proposed,3));
% sys0 = Dynamics.LotkaVolterraContNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,(k-1)*2+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h4 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,(k-1)*2+1) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

load('cont_100.mat')
idx = 1:20:size(Q_proposed,3);
% col = cmap(5,:);
col = 'c';
% col = Utils.getColor(5, [1,5], cmap);
% t = linspace(0,1,size(Q_proposed,3));
% sys0 = Dynamics.LotkaVolterraContNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
for k = idx
    if k > 1
        %%% Proposed method
        x_ = sys0.xN(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts).x;
        h5 = plot(x_(1,:), x_(2,:), 'o', 'color', col, 'linewidth', 1);
        x_ = sys0.xN(:,k) + Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
%         x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,nPts2).x;
        plot(x_(1,:), x_(2,:), '-', 'color', col, 'linewidth', 1);
    end
end

k = 1;
text(sys0.xN(1,k), sys0.xN(2,k), ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 21;
text(sys0.xN(1,k), sys0.xN(2,k)-0.13, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 41;
text(sys0.xN(1,k), sys0.xN(2,k)-0.17, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 61;
text(sys0.xN(1,k), sys0.xN(2,k)-0.16, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 81;
text(sys0.xN(1,k), sys0.xN(2,k)-0.19, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)
k = 101;
text(sys0.xN(1,k), sys0.xN(2,k)-0.21, ['$t=',num2str(t(k)),'$ s'], 'horizontalalignment', 'center', 'fontsize', 14)

legend([h0,h1,h2,h3,h4,h5],...
    'HJB PDE',...
    '$N = 5$',...
    '$N = 10$',...
    '$N = 20$',...
    '$N = 50$',...
    '$N = 100$',...
    'location', 'northeast')

xlabel('$x_1$')
ylabel('$x_2$')
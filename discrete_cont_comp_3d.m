clear; close all; clc;

%% 3D visualization
t = linspace(0,1,101);
figure;
cla; hold on; grid on; axis tight;
set(gcf, 'position', [681,559,480,420])
h7 = Utils.plot_set(X, t, length(t), 'k', 0.3);
% h9 = Utils.plot_set(F_sos, t, length(t), 'r', 0.3);
% h8 = Utils.plot_set(F_proposed, t, length(t), [65,169,76]/255, 0.3);

for k = round(linspace(1,length(t),5))
    if k > 1
        %%% SOS
        x_ = Q_sos(:,:,k)^(1/2)*Math.Sphere(1,500).x;
        h3 = plot3(x_(1,:), t(k)*ones(size(x_(1,:))), x_(2,:), '-', 'color', 'r', 'linewidth', 2);
        %%% Proposed
        x_ = Q_proposed(:,:,k)^(1/2) * Math.Sphere(1,500).x;
        h5 = plot3(x_(1,:), t(k)*ones(size(x_(1,:))), x_(2,:), '-', 'color', [65,169,76]/255, 'linewidth', 2);
    end
    %%% HJB eqn
    x_ = X{k};
    h6 = plot3(x_(1,:), t(k)*ones(size(x_(1,:))), x_(2,:), 'k-', 'linewidth', 2);
end
xlabel('$\tilde{x}_1$')
ylabel('$t$ [s]')
zlabel('$\tilde{x}_2$')
view([122,11])
xlim([-0.1,0.1])
zlim([-0.1,0.1])
ylim([t(1), t(end)])
set(gca, 'xdir', 'reverse')
camlight('left')
legend([h7,h9,h8],...
    'HJB PDE',...
    'SOS program',...
    'Proposed',...
    'location', 'northwest')
camlight('right')
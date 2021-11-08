% clear;
close all; clc;

addpath(genpath('../reachability_analysis/3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('../reachability_analysis/3rd_party/ToolboxLS'))

n = 2;
xlim = [0.75, 1.3];
ylim = [0.85, 1.25];

K = [7, 7]; % number of regions
xgr = cell(2,1);
xgr{1} = linspace(xlim(1), xlim(end), K(1)+1);
xgr{2} = linspace(ylim(1), ylim(end), K(2)+1);

Null = cell(2,1);
A = cell(2,1);
for i = 1:(length(xgr{1})-2)
%     A_ = [1, xgr{1}(i+1), 0; 0, 0, 1];
%     A_ = [1, xgr{1}(i+1), 0];
    A_ = [0, 0, 1];
    A{1}(:,:,i) = A_;
    Null{1}(:,:,i) = null(A_);
end
for i = 1:(length(xgr{2})-2)
%     A_ = [1, 0, xgr{2}(i+1); 0, 1, 0];
%     A_ = [1, 0, xgr{2}(i+1)];
    A_ = [0, 1, 0];
    A{2}(:,:,i) = A_;
    Null{2}(:,:,i) = null(A_);
end

I_region = Polynomial.integer_grid(K-1);
id_region = Polynomial.getId(I_region, K-1);

t = linspace(0,0.2,101);

%% Target function
gr_all = createGrid([xlim(1),ylim(1)],[xlim(2),ylim(2)],[51,51]);
Ntarget = [3, 3, 4];
Itarget = Polynomial.integer_grid(Ntarget);
htarget = prod((1./(Itarget+1)).*( ones(size(Itarget)).^(Itarget+1) - (-ones(size(Itarget))).^(Itarget+1) ),2);
% ctarget = randn(size(Itarget,1),1);
ptarget = Polynomial(ctarget,Itarget);
Btarget = Polynomial.bernstein_transform_matrix(Ntarget);

%% Approximation function
N = [1,1,3];
I = Polynomial.integer_grid(N);
effective = sum(I(:,1:2),2) < 2;
I_effective = I( effective, : );
id_effective = Polynomial.getId( I_effective, Ntarget );

proj = zeros(size(Itarget,1), length(id_effective));
proj(id_effective,:) = eye(length(id_effective));

%% Matrices per grid
TFtarget = cell( prod(K), 1 );
gr_grid = cell( prod(K), 1 );
for i = 1:length(id_region)
    id_ = id_region(i);
    order_ = Polynomial.getOrder( id_, K-1 );
    lb_ = [xgr{1}(order_(1)+1); xgr{2}(order_(2)+1); t(1)];
    ub_ = [xgr{1}(order_(1)+2); xgr{2}(order_(2)+2); t(end)];
    TFtarget{id_} = Polynomial.domain_transform_matrix(Ntarget, (ub_-lb_)*0.5, (ub_+lb_)*0.5);
    gr_grid{id_} = createGrid(lb_(1:2), ub_(1:2), [2,2]);
end

%% Optimization per grid
cvx_begin
    variable c( size(I_effective,1), length(id_region) )
    variable v1( size(Null{1},2)*(N(end)+1), K(1)-1 )
    variable v2( size(Null{1},2)*(N(end)+1), K(2)-1 )
    
    cost = 0.0;
    for k = 1:length(id_region)
        cost = cost + -htarget'*TFtarget{k}*proj*c(:,k);
        (Btarget\TFtarget{k})*(ctarget - proj*c(:,k)) >= 0
    end
    
    visited = zeros(length(id_region), 1);
    visited( id_region(1) ) = true;
    leaf = id_region(1);
    
    iter = 1;
    while (~prod(visited))% && iter < 3
        new_leaf_ = [];
        for i = 1:length(leaf)
            id_ = leaf(i);
            order_ = Polynomial.getOrder( id_, K-1 );
    
            for j = 1:n
                increment_ = zeros(1,n);
                increment_(j) = 1.0;
                order_adj_ = order_ + increment_;
                id_adj_ = Polynomial.getId( order_adj_, K-1 );
                
                if id_adj_ <= id_region(end)
                    if ~visited(id_adj_) 
                        v_id_ = order_adj_(j);
                        if j == 1
%                             reshape(c(:,id_adj_),[3,N(end)+1]) == reshape(c(:,id_),[3,N(end)+1]) +...
%                                 Null{j}(:,:,order_adj_(j))*reshape(v1(:,v_id_),[size(Null{j}(:,:,order_adj_(j)),2),N(end)+1]);
                        elseif j == 2
%                             reshape(c(:,id_adj_),[3,N(end)+1]) == reshape(c(:,id_),[3,N(end)+1]) +...
%                                 Null{j}(:,:,order_adj_(j))*reshape(v2(:,v_id_),[size(Null{j}(:,:,order_adj_(j)),2),N(end)+1]);
                        end
                    
                        new_leaf_ = cat(1, new_leaf_, id_adj_);
                        visited(id_adj_) = true;
                    end
                end
            end
        end
        leaf = new_leaf_;
        iter = iter + 1;
    end

    minimize cost
cvx_end
    
%%
figure(1)
hold on; grid on; axis tight;
xlabel('$x_1$')
ylabel('$x_2$')
ax = gca;
ax.XLim = xlim;
ax.YLim = ylim;
for i = round(linspace(1,length(t),11))%1:length(t)%
    cla;
    t_ = t(i);
    t_basis_ = [1; t_; t_^2; t_^3];
    for k = 1:length(id_region)
        C_ = reshape(c(:,k), [3,N(end)+1]);
        coeff_ = C_ * t_basis_;
        V_ = coeff_(1) + gr_grid{k}.xs{1}*coeff_(2) + gr_grid{k}.xs{2}*coeff_(3);
        surf(gr_grid{k}.xs{1}, gr_grid{k}.xs{2}, V_, 'facealpha', 0.5, 'facecolor', 'r')
    end
    
    p_ = Polynomial(ctarget.*t_.^ptarget.order(:,end), ptarget.order(:,1:end-1));
    V_ = p_.eval(gr_all);
    surf(gr_all.xs{1}, gr_all.xs{2}, V_, 'linestyle', 'none', 'facealpha', 0.5)
    view([31,15])
    drawnow
end


%% Target validation


% Nt = 4; % time order

% c0 = rand(3,Nt+1);
% v = cell(2,1);
% v{1} = randn(length(xgr{1})-2,Nt+1);
% v{2} = randn(length(xgr{2})-2,Nt+1);

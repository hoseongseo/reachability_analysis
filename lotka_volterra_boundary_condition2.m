clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))

wMax = 0.1;
%%% system (polynomial dynamics)
t = linspace(0,1,51);
sys0 = Dynamics.LotkaVolterraNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
sys = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin
% sys = Dynamics.LotkaVolterra2( wMax );

% polynomial dynamics
xp(sys.Nx,1) = Polynomial;
for i = 1:sys.Nx
    ord_ = zeros(1,sys.Nx+sys.Nw);
    ord_(i) = 1;
    xp(i) = Polynomial(1, ord_);
end
wp(sys.Nw,1) = Polynomial;
for i = 1:sys.Nw
    ord_ = zeros(1,sys.Nx+sys.Nw);
    ord_(sys.Nx+i) = 1;
    wp(i) = Polynomial(1, ord_);
end

fp = sys.f(xp, zeros(sys.Nu,1), wp, t(1));
fp_order_all = [];
for i = 1:sys.Nx
    fp_order_all = [fp_order_all; fp(i).order];
end
N_f = max(reshape(fp_order_all, [numel(fp_order_all)/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);


% approximated value function
N = [2, 2]; % order = x1, x2
% N_a = [N(1:sys.Nx), zeros(1,sys.Nw), N(end)]; % order = x1, x2, w, t
N_a = [N, zeros(1,sys.Nw)]; % augmentation
N_max = N_a + N_f;
B_N_max = Polynomial.bernstein_transform_matrix(N_max);
P_N_N_f = Polynomial.product_matrix(N_a, N_f);
D = Polynomial.derivative_matrix(N_a);
M = zeros( prod(N_a+N_f+1), prod(N+1), length(t)-1 );
for k = 1:length(t)-1
    mat = zeros(prod(N_a+N_f+1), prod(N+1));
    f_ = sys.f(xp, zeros(sys.Nu,1), wp, t(k));
    for j = 1:sys.Nx
        tmp = P_N_N_f * kron( eye(prod(N+1)), Polynomial.expand_matrix(N_f, f_(j).order)*f_(j).coeff ) * D(:,:,j);
        mat = mat + tmp;
    end
    M(:,:,k) = mat;
end

E_N = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a));
E_N_N_f = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a + N_f));
I_N_a = Polynomial.integer_grid( N_a );
h_N_a = prod((1./(I_N_a+1)).*( ones(size(I_N_a)).^(I_N_a+1) - (-ones(size(I_N_a))).^(I_N_a+1) ),2);


% initial coefficient
q0 = [1.2; 1.1];
Q0 = diag([0.05; 0.05])^2;
% V0 = (xp-q0)'*inv(Q0)*(xp-q0) - 1;
V0 = (xp)'*inv(Q0)*(xp) - 1;
c0 = Polynomial.expand_matrix(N(1:sys.Nx), V0.order(:,1:sys.Nx))*V0.coeff;

%%
% region subdivision
% lb = [0.6; 0.7];
% ub = [1.3; 1.3];
lb = [-0.08; -0.08];
ub = [0.08; 0.08];
K = [3, 3];

% Grid
xgr = cell(1,length(N));
for j = 1:length(N)
    xgr{j} = linspace(lb(j), ub(j), K(j)+1);
end

% Number of additional variables
nvar = zeros(1,length(N));
for j = 1:length(N)
    jc = setdiff(1:length(N), j);
    
%     if (j <= sys.Nx + sys.Nw) && (j >= sys.Nx + 1)
%         nvar(j) = 0; % for disturbance, continuity must not be considered
%     else
        nvar(j) = prod(N+1) - 2*prod(N(jc)+1); % total - continuous constraint
%     end
end

% Nullspace gain along grid
n = length(N);
Null = cell(1,n);
I = Polynomial.integer_grid( N );
for j = 1:n
    Null{j} = zeros(prod(N+1), nvar(j), K(j)-1);
    
    jc = setdiff(1:n, j);
    Nc = N(jc);
    Ic = Polynomial.integer_grid(Nc);
    idc = Polynomial.getId( I(:,jc), Nc );
    
%     if (j <= sys.Nx + sys.Nw) && (j >= sys.Nx + 1)
%         % do nothing
%     else
        for k = 1:(K(j)-1)
            C_ = zeros(prod(N+1), prod(Nc+1));
            for i = 1:size(C_,1)
                C_(i,idc(i)) = prod( xgr{j}(k+1).^I(i,j));
            end
            L_ = null( [C_, D(:,:,j)'*C_]' );
            Null{j}(:,:,k) = L_;
        end
%     end
end

% domain transform matrix
I_region = Polynomial.integer_grid(K-1);
id_region = Polynomial.getId(I_region, K-1);
T_max = zeros( prod(N_max+1), prod(N_max+1), length(id_region) );
T_a = zeros( prod(N_a+1), prod(N_a+1), length(id_region) );
for k = 1:length(id_region)
    K_ = I_region(k,:);
    
    lb_ = zeros(n,1);
    ub_ = zeros(n,1);
    for j = 1:n
        lb_(j) = xgr{j}(K_(j)+1);
        ub_(j) = xgr{j}(K_(j)+2);
    end
    a_ = 0.5*(ub_-lb_);
    b_ = 0.5*(ub_+lb_);
    
    T_max(:,:,k) = Polynomial.domain_transform_matrix(N_max, [a_; 1], [b_; 0]);
    T_a(:,:,k) = Polynomial.domain_transform_matrix(N_a, [a_; 1], [b_; 0]);
end

c_prev = repmat(c0, [1, prod(K)]);

%% Node selection

% for k = 1:(length(t)-1)

cvx_begin
variable c_test( prod(N+1), prod(K), length(t) )
variable v_test( max(nvar), sum(K-1), length(t) )

c_test(:,:,1) == c_prev
v_test(:,:,1) == zeros(max(nvar), sum(K-1))
cost = 0.0;

for k = 2:(length(t))

visited = zeros(length(id_region), 1);
visited( id_region(1) ) = true;
leaf = id_region(1);

(B_N_max \ T_max(:,:,1)) * ( E_N*(c_test(:,1,k) - c_test(:,1,k-1)) + (t(k)-t(k-1))*(E_N_N_f*M(:,:,k-1)*c_test(:,1,k-1)) ) <= 0 % constraint for the initial region
% (B_N_max \ T_max(:,:,1)) * (E_N * D(:,:,sys.Nx+1) + M) * c_test(:,1) <= 0 % constraint for the initial region
% c_test(1:prod(N(1:sys.Nx)+1),1) == c0; % IC per region
cost = cost - (1/prod(K))*h_N_a'*(T_a(:,:,1)*c_test(:,1,k)); % cost per region

while ~prod(visited)
    new_leaf_ = [];
    for i = 1:length(leaf)
        id_ = leaf(i);
        order_ = Polynomial.getOrder( id_, K-1 );
        
        for j = 1:n
            increment_ = zeros(1,n);
            increment_(j) = 1.0;
            
            order_adj_ = order_ + increment_;
            id_adj_ = Polynomial.getId( order_adj_, K-1 );
            
            if id_adj_ <= id_region(end) % consider the maximum id
                if visited(id_adj_)
                    % DO NOTHING
                    % disp("Already visited")
                else
                    % DO SOMETHING
                    v_order_ = [j, order_adj_(j)];
                    v_id_ = order_adj_(j);
                    if j > 1
                        v_id_ = v_id_ + sum( K(1:(j-1))-1 );
                    end
%                     c_test(:,id_adj_) = c_test(:,id_) + Null{j}(:,:,order_adj_(j)) * v_test(:,v_id_);                    
                    % disp( ['c_',num2str(id_adj_),' is computed from c_',num2str(id_),' with v(',num2str(j),',',num2str(order_adj_(j)),') and v_id = ',num2str(v_id_)] )
                    
                    % apply constraints
                    c_test(:,id_adj_,k) == c_test(:,id_,k) + Null{j}(:,:,order_adj_(j)) * v_test(1:nvar(j),v_id_,k)
%                     B_N_max \ T_max(:,:,id_adj_) * (E_N * D(:,:,sys.Nx+1) + M) * c_test(:,id_adj_) <= 0
%                     c_test(1:prod(N(1:sys.Nx)+1),id_adj_) == c0
                    (B_N_max \ T_max(:,:,id_adj_)) * ( E_N*(c_test(:,id_adj_,k) - c_test(:,id_adj_,k-1)) + (t(k)-t(k-1))*(E_N_N_f*M(:,:,k-1)*c_test(:,id_adj_,k-1)) ) <= 0
                    % update cost
                    cost = cost - (1/prod(K))*h_N_a'*(T_a(:,:,id_adj_)*c_test(:,id_adj_,k)); % cost per region
                    
                    new_leaf_ = cat(1, new_leaf_, id_adj_); % update the next leaf
                    visited(id_adj_) = true; % trun on flag such that the corresponding region is visited
                end
            end
        end
    end
    leaf = new_leaf_;
end

end

minimize cost
cvx_end


%% visualization
X_hist = cell(1,length(t));
X_hist{1} = Q0^(1/2) * Math.Sphere(1,200).x;

for k = 2:length(t)
figure(121)
hold on; grid on; axis tight; 
cla
X0 = [];
for kk = 1:prod(K)
    K_ = I_region(kk,:);
    
    lb_ = zeros(n,1);
    ub_ = zeros(n,1);
    for j = 1:n
        lb_(j) = xgr{j}(K_(j)+1);
        ub_(j) = xgr{j}(K_(j)+2);
    end
    gr_ = createGrid( lb_, ub_, 101*ones(size(lb_)) );
    p_ = Polynomial( c_test(:,kk,k), I );
    V_ = p_.eval( gr_ );
    X_ = Utils.get_level_set(gr_, V_, 0.0);
    X0 = [X0, X_];
    surf( gr_.xs{1}, gr_.xs{2}, V_, 'linestyle', 'none', 'facealpha', 0.5, 'facecolor', rand(3,1))
    
    surf(gr.xs{1}, gr.xs{2}, V(:,:,k), 'linestyle', 'none', 'facealpha', 0.5)
end
plot3(X0(1,:), X0(2,:), zeros(1,size(X0,2)), 'k.')
X_hist{k} = X0;
% view([0,90])
view([17,16])
xlabel('$x_1$')
ylabel('$x_2$')
xlim([lb(1),ub(1)])
ylim([lb(2),ub(2)])
drawnow
end

% update
% c_prev = c_test;
% end

%% HJB
nGrid = [101; 101];
minGrid = lb;
maxGrid = ub;
gr = createGrid(minGrid, maxGrid, nGrid);

V0 = gr.xs{1}.*gr.xs{1}/Q0(1,1) + gr.xs{2}.*gr.xs{2}/Q0(2,2) - 1;
X0 = Utils.get_level_set(gr, V0, 0.0);

% solve
hjb_equation = HJBequation(sys, gr);
V = zeros([size(V0), length(t)]);
V(:,:,1) = V0;
tic
for i = 1:length(t)-1
    V(:,:,i+1) = hjb_equation.solve(V(:,:,i), t(i), t(i+1));
end
toc

% extract zero-level set
X = cell(1,length(t));
X{1} = X0;
for i = 2:length(t)
    X{i} = Utils.get_level_set(gr, V(:,:,i), 0.0);
end

%%
figure(1233)
cla; hold on; grid on; axis tight; 
for k = round(linspace(1,length(t),21))
    tmp = sys0.xN(:,k) + X{k};
    plot(tmp(1,:), tmp(2,:), 'k-', 'linewidth', 2)
    
    tmp = sys.xN(:,k) + X_hist{k};
    plot(tmp(1,:), tmp(2,:), 'b.', 'linewidth', 2)
end
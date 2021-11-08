clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))

wMax = 0.1;
sys = Dynamics.LotkaVolterra2( wMax );

% polynomial dynamics
xp(sys.Nx,1) = Polynomial;
for i = 1:sys.Nx
    ord_ = zeros(1,sys.Nx+sys.Nw+1);
    ord_(i) = 1;
    xp(i) = Polynomial(1, ord_);
end
wp(sys.Nw,1) = Polynomial;
for i = 1:sys.Nw
    ord_ = zeros(1,sys.Nx+sys.Nw+1);
    ord_(sys.Nx+1+i) = 1;
    wp(i) = Polynomial(1, ord_);
end
tp = Polynomial(1, [zeros(1,sys.Nx),1,zeros(1,sys.Nw)]); % x1, x2, t, w

fp = sys.f(xp, zeros(sys.Nu,1), wp, tp);
fp_order_all = [];
for i = 1:sys.Nx
    fp_order_all = [fp_order_all; fp(i).order];
end
N_f = max(reshape(fp_order_all, [numel(fp_order_all)/(sys.Nx+sys.Nw+1), (sys.Nx+sys.Nw+1)]), [], 1);

% approximated value function
N = [3, 3, 11]; % order = x1, x2, t
% N_a = [N(1:sys.Nx), zeros(1,sys.Nw), N(end)]; % order = x1, x2, w, t
N_a = [N, 0]; % augmented order = x1, x2, t, w
N_max = N_a + N_f;
B_N_max = Polynomial.bernstein_transform_matrix(N_max);
P_N_N_f = Polynomial.product_matrix(N_a, N_f);
D = Polynomial.derivative_matrix(N_a);
M = zeros(prod(N_max+1), prod(N_a+1));
for j = 1:sys.Nx
    M = M + P_N_N_f * kron( eye(prod(N_a+1)), Polynomial.expand_matrix(N_f, fp(j).order)*fp(j).coeff ) * D(:,:,j);
end
E_N = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a));
I_N_a = Polynomial.integer_grid( N_a );
h_N_a = prod((1./(I_N_a+1)).*( ones(size(I_N_a)).^(I_N_a+1) - (-ones(size(I_N_a))).^(I_N_a+1) ),2);

% initial coefficient
q0 = [1.2; 1.1];
Q0 = diag([0.05; 0.05])^2;
V0 = (xp-q0)'*inv(Q0)*(xp-q0) - 1;
c0 = Polynomial.expand_matrix(N(1:sys.Nx), V0.order(:,1:sys.Nx))*V0.coeff;


% region subdivision
% lb = [1.0; 1.0; 0]; % x1, x2, t
% ub = [1.3; 1.3; 0.25];
lb = [0.6; 0.7; 0]; % x1, x2, t
ub = [1.3; 1.3; 1.0];
K = [2, 2, 2];
% n = length(N);

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


%% Node selection

cvx_begin
variable c_test( prod(N+1), prod(K) )
variable v_test( max(nvar), sum(K-1) )

visited = zeros(length(id_region), 1);
visited( id_region(1) ) = true;
leaf = id_region(1);

(B_N_max \ T_max(:,:,1)) * (E_N * D(:,:,sys.Nx+1) + M) * c_test(:,1) <= 0 % constraint for the initial region
c_test(1:prod(N(1:sys.Nx)+1),1) == c0; % IC per region
cost = - h_N_a'*(T_a(:,:,1)*c_test(:,1)); % cost per region

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
                    c_test(:,id_adj_) == c_test(:,id_) + Null{j}(:,:,order_adj_(j)) * v_test(1:nvar(j),v_id_)
                    B_N_max \ T_max(:,:,id_adj_) * (E_N * D(:,:,sys.Nx+1) + M) * c_test(:,id_adj_) <= 0
                    c_test(1:prod(N(1:sys.Nx)+1),id_adj_) == c0
                    
                    % update cost
                    cost = cost - h_N_a'*(T_a(:,:,id_adj_)*c_test(:,id_adj_)); % cost per region
                    
                    new_leaf_ = cat(1, new_leaf_, id_adj_); % update the next leaf
                    visited(id_adj_) = true; % trun on flag such that the corresponding region is visited
                end
            end
        end
    end
    leaf = new_leaf_;
end

minimize cost
cvx_end


%% visualization
jc = setdiff(1:n, n);

Nc = N(jc);
Ic = Polynomial.integer_grid(Nc);
idc = Polynomial.getId( I(:,jc), Nc );
C = zeros(prod(N+1), prod(Nc+1));


% Ntargetc = Ntarget(jc);
% Itargetc = Polynomial.integer_grid(Ntargetc);
% idtargetc = Polynomial.getId( Itarget(:,jc), Ntargetc );
% Ctarget = zeros(prod(Ntarget+1), prod(Ntargetc+1));
% 
% grtarget = createGrid([-1,-1],[1,1],[51,51]);

figure(121)
hold on; grid on; axis tight; 

for t = linspace(lb(end), ub(end), 101)
% for t = 0
cla;

% approximated function
for i = 1:size(C,1)
    C(i,idc(i)) = prod( t.^I(i,end)); % C' * c_test = polynomial of state at the quired time
end
% get region idx according to time
t_idx = 1;
for i = 2:length(xgr{end})
    if t > xgr{end}(i)
        t_idx = t_idx + 1;
    end    
end
region_idx_ = id_region( prod(K(1:end-1))*(t_idx-1) + (1:prod(K(1:end-1))) );
X0 = [];
for k = 1:prod(K(1:end-1))
    K_ = I_region(k,:);
    
    lb_ = zeros(n-1,1);
    ub_ = zeros(n-1,1);
    for j = 1:n-1
        lb_(j) = xgr{j}(K_(j)+1);
        ub_(j) = xgr{j}(K_(j)+2);
    end
    gr_ = createGrid( lb_, ub_, 101*ones(size(lb_)) );
    p_ = Polynomial( C'*c_test(:,region_idx_(k)), Ic );
    V_ = p_.eval( gr_ );
    X_ = Utils.get_level_set(gr_, V_, 0.0);
    X0 = [X0, X_];
    surf( gr_.xs{1}, gr_.xs{2}, V_, 'linestyle', 'none', 'facealpha', 0.5, 'facecolor', rand(3,1))
end
plot3(X0(1,:), X0(2,:), zeros(1,size(X0,2)), 'k.')

% % target function
% for i = 1:size(Ctarget,1)
%     Ctarget(i,idtargetc(i)) = prod( t.^Itarget(i,end)); % C' * c_test = polynomial of state at the quired time
% end
% ctarget_ = Ctarget' * ctarget;
% ptarget_ = Polynomial(ctarget_, Itargetc);
% surf(grtarget.xs{1}, grtarget.xs{2}, ptarget_.eval( grtarget ), 'linestyle', '-', 'facealpha', 0.2, 'facecolor', 'none')

view([0,90])
xlabel('$x_1$')
ylabel('$x_2$')
% pause
drawnow
end

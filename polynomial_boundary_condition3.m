clear; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))

N = [1, 1, 3];
n = length(N);
I = Polynomial.integer_grid(N);
h = prod((1./(I+1)).*( ones(size(I)).^(I+1) - (-ones(size(I))).^(I+1) ),2);
D = Polynomial.derivative_matrix(N);

proj = zeros(prod(N+1),3*4);
proj([1,2,3,5,6,7,9,10,11,13,14,15],:) = eye(12);


lb = [-1; -1; 0];
ub = [1; 1; 2];

K = [11, 11, 1];

% Grid
xgr = cell(1,n);
for j = 1:n
    xgr{j} = linspace(lb(j), ub(j), K(j)+1);
end

% Number of additional variables
nvar = zeros(1,n);
for j = 1:n
    jc = setdiff(1:n, j);
    nvar(j) = prod(N+1) - prod(N(jc)+1); % total - continuous constraint
end

% Nullspace gain along grid
Null = cell(1,n);
for j = 1:n
    Null{j} = zeros(prod(N+1), nvar(j), K(j)-1);
    
    jc = setdiff(1:n, j);
    Nc = N(jc);
    Ic = Polynomial.integer_grid(Nc);
    idc = Polynomial.getId( I(:,jc), Nc );
    for k = 1:(K(j)-1)
        C_ = zeros(prod(N+1), prod(Nc+1));
        for i = 1:size(C_,1)
            C_(i,idc(i)) = prod( xgr{j}(k+1).^I(i,j));
        end
%         L_ = null( [C_, D(:,:,j)'*C_]' );
        L_ = null( C_' );
        Null{j}(:,:,k) = L_;
    end
end

I_region = Polynomial.integer_grid(K-1);
id_region = Polynomial.getId(I_region, K-1);

%%
% target function
Ntarget = [3,3,3];
Itarget = Polynomial.integer_grid(Ntarget);
Btarget = Polynomial.bernstein_transform_matrix(Ntarget);
Etarget = Polynomial.expand_matrix(Ntarget, I);
htarget = prod((1./(Itarget+1)).*( ones(size(Itarget)).^(Itarget+1) - (-ones(size(Itarget))).^(Itarget+1) ),2);

% initial coefficient
ctarget = randn(prod(Ntarget+1), 1);
% x10 = -0.2;
% x20 = 0.3;
% r0 = 0.5;
% ctarget = zeros(prod(Ntarget+1), 1);
% ctarget(1) = x10^2 + x20^2 - r0^2;
% ctarget(2) = -2*x10;
% ctarget(6) = -2*x20;
% ctarget(3) = 1.0;
% ctarget(11) = 1.0;
% ctarget(4) = -1.1637;
% ctarget(5) = randn(1);


% Domain setting
T = zeros( prod(Ntarget+1), prod(Ntarget+1), length(id_region) );
% gr = cell(1,length(id_region));
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
    
    T(:,:,k) = Polynomial.domain_transform_matrix(Ntarget, a_, b_);
%     gr{k} = createGrid( lb_, ub_, 101*ones(size(lb_)) );
end

%% Node selection

% c_test = zeros( prod(N+1), prod(K) );
% c_test(:,1) = randn( prod(N+1),1 );
% v_test = randn(nvar, sum(K-1));

cvx_begin
% variable c_test( prod(N+1), prod(K) )
variable c_test( size(proj,2), prod(K) )
variable v_test( max(nvar), sum(K-1) )

% c_test = zeros(prod(N+1), prod(K));
% c_test(:,1) = c0;

visited = zeros(length(id_region),1);
visited( id_region(1) ) = true;
leaf = id_region(1);

% constraint and cost on initial domain
% Btarget \ T(:,:,1) * (ctarget - Etarget * c_test(:,1)) >= 0
% cost = - htarget'*(T(:,:,1)*Etarget*c_test(:,1));
Btarget \ T(:,:,1) * (ctarget - Etarget * proj * c_test(:,1)) >= 0
cost = - htarget'*(T(:,:,1)*Etarget * proj * c_test(:,1));

% node iteration
while ~prod(visited) % repeat this until all regions are visited
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
%                     v_order_ = [j, order_adj_(j)];
                    v_id_ = order_adj_(j);
                    if j > 1
                        v_id_ = v_id_ + sum( K(1:(j-1))-1 );
                    end
%                     c_test(:,id_adj_) = c_test(:,id_) + Null{j}(:,:,order_adj_(j)) * v_test(:,v_id_);                    
                    % disp( ['c_',num2str(id_adj_),' is computed from c_',num2str(id_),' with v(',num2str(j),',',num2str(order_adj_(j)),') and v_id = ',num2str(v_id_)] )
                    
%                     % apply constraints
%                     c_test(:,id_adj_) == c_test(:,id_) + Null{j}(:,:,order_adj_(j)) * v_test(1:nvar(j),v_id_)
%                     Btarget \ T(:,:,id_adj_) * (ctarget - Etarget * c_test(:,id_adj_)) >= 0
%                     
%                     % update cost
%                     cost = cost - htarget'*(T(:,:,id_adj_)*Etarget*c_test(:,id_adj_));
                    
                    % apply constraints
                    proj*c_test(:,id_adj_) == proj*c_test(:,id_) + Null{j}(:,:,order_adj_(j)) * v_test(1:nvar(j),v_id_)
                    Btarget \ T(:,:,id_adj_) * (ctarget - Etarget * proj * c_test(:,id_adj_)) >= 0
                    
                    % update cost
                    cost = cost - htarget'*(T(:,:,id_adj_)*Etarget*proj*c_test(:,id_adj_));
                    
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


Ntargetc = Ntarget(jc);
Itargetc = Polynomial.integer_grid(Ntargetc);
idtargetc = Polynomial.getId( Itarget(:,jc), Ntargetc );
Ctarget = zeros(prod(Ntarget+1), prod(Ntargetc+1));

grtarget = createGrid([-1,-1],[1,1],[51,51]);



figure(121)
hold on; grid on;

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
for k = 1:prod(K(1:end-1))
    K_ = I_region(k,:);
    
    lb_ = zeros(n-1,1);
    ub_ = zeros(n-1,1);
    for j = 1:n-1
        lb_(j) = xgr{j}(K_(j)+1);
        ub_(j) = xgr{j}(K_(j)+2);
    end
    gr_ = createGrid( lb_, ub_, 2*ones(size(lb_)) );
%     p_ = Polynomial( C'*c_test(:,region_idx_(k)), Ic );
    p_ = Polynomial( C'*proj*c_test(:,region_idx_(k)), Ic );
    surf( gr_.xs{1}, gr_.xs{2}, p_.eval( gr_ ), 'linestyle', 'none', 'facealpha', 0.5, 'facecolor', rand(3,1))
end

% target function
for i = 1:size(Ctarget,1)
    Ctarget(i,idtargetc(i)) = prod( t.^Itarget(i,end)); % C' * c_test = polynomial of state at the quired time
end
ctarget_ = Ctarget' * ctarget;
ptarget_ = Polynomial(ctarget_, Itargetc);
surf(grtarget.xs{1}, grtarget.xs{2}, ptarget_.eval( grtarget ), 'linestyle', '-', 'facealpha', 0.2, 'facecolor', 'none')

view([-29,20])
xlabel('$x_1$')
ylabel('$x_2$')
% pause
drawnow
end

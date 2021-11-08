clear; close all; clc;

% n = 3;
% N = [2,2,2];
% I = Polynomial.integer_grid(N);
% j = 1;
% xj = 0.5;
% D = Polynomial.derivative_matrix(N);
% 
% jc = setdiff(1:n, j);
% Nc = N(jc);
% Ic = Polynomial.integer_grid(Nc);
% 
% idc = Polynomial.getId( I(:,jc), Nc );
% C = zeros(prod(N+1), prod(Nc+1));
% for i = 1:size(C,1)
%     C(i,idc(i)) = prod(xj.^I(i,j));
% end
% L = null([C, D(:,:,j)'*C]');

%% region spacification
N = [2, 2, 2];
n = length(N);
I = Polynomial.integer_grid(N);
D = Polynomial.derivative_matrix(N);

lb = [-0.5; -0.5; -0.5]; % lower bound
ub = [1.5; 1.5; 1.5];   % upper bound
K = [4, 3, 2];        % sub-division parameter

% Grid
xgr = cell(1,n);
for j = 1:n
    xgr{j} = linspace(lb(j), ub(j), K(j)+1);
end

% Number of additional variables
nvar = zeros(1,n);
for j = 1:n
    jc = setdiff(1:n, j); % state index (complement)
    Nc = N(jc); % degree (complement) 
    nvar(j) = (prod(N+1) - 2*prod(Nc+1));
end

% Nullspace gain along grid
Nullgr = cell(1,n);
for j = 1:n
    
    jc = setdiff(1:n, j); % state index (complement)
    Nc = N(jc); % degree (complement)
    
%     nvar = (prod(N+1) - 2*prod(Nc+1))*(K(j)-1);
%     Nullgr{j} = zeros( prod(N+1), prod(N+1) - 2*prod(Nc+1), K(j)-1 ); % nullspace gain per each state axis
    Nullgr{j} = cell(1,K(j)-1);

%     Nullgr{j} = [];
    
    Ic = Polynomial.integer_grid(Nc);
    idc = Polynomial.getId( I(:,jc), Nc );
    
    for k = 1:(K(j)-1)
        C_ = zeros(prod(N+1), prod(Nc+1));
        for i = 1:size(C_,1)
            C_(i,idc(i)) = prod( xgr{j}(k+1).^I(i,j));
        end
        L_ = null( [C_, D(:,:,j)'*C_]' );
%         Nullgr{j} = [Nullgr{j}, L_];
%         Nullgr{j}(:,:,k) = L_;
        Nullgr{j}{k} = L_;
    end
end

%%
% Regions
I_region = Polynomial.integer_grid(K-1);
id_region = Polynomial.getId(I_region, K-1);


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
    
    
end

% Nullspace gain along idx
% Null_idx = zeros( prod(N+1), sum( nvar.*(K-1) ), length(id_region) );
% for i = 1:length(id_region)
%     I_region_ = I_region(i,:);
%     for j = 1:n
%         for k = 
%         
%         end
%     end
% end



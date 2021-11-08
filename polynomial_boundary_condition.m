clear; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))

N = [2,2];
n = length(N);
I = Polynomial.integer_grid(N);
h = prod((1./(I+1)).*( ones(size(I)).^(I+1) - (-ones(size(I))).^(I+1) ),2);
D = Polynomial.derivative_matrix(N);

Du = Utils.duplication_matrix(n);

Cslice1 = zeros(N(2)+1, size(I,1)); % x1 = const 
for i = 1:size(Cslice1,1)
    Cslice1(i,:) = (I(:,2) == (i-1));
end
Cslice2 = zeros(N(1)+1, size(I,1)); % x2 = const
for i = 1:size(Cslice2,1)
    Cslice2(i,:) = (I(:,1) == (i-1));
end

% in this simulation, the whole range is [-1,1]
x1gr = linspace(-1,1,7);
x2gr = linspace(-1,1,7);

x1Nu = cell( length(x1gr)-2, 1 );
for i = 1:size(x1Nu,1)
    C_ = Cslice1 * diag( x1gr(i+1).^I(:,1) );
    dC_ = C_ * D(:,:,1);
    Ctotal_ = [C_; dC_];
    x1Nu{i} = null(Ctotal_);
%     x1Nu{i} = null(C_);
end

x2Nu = cell( length(x2gr)-2, 1 );
for i = 1:size(x2Nu,1)
    C_ = Cslice2 * diag( x2gr(i+1).^I(:,2) );
    dC_ = C_ * D(:,:,2);
    Ctotal_ = [C_; dC_];
    x2Nu{i} = null(Ctotal_);
%     x2Nu{i} = null(C_);
end

Ntarget = [4,4];
Itarget = Polynomial.integer_grid(Ntarget);
% ctarget = randn(prod(Ntarget+1), 1);
x10 = -0.2;
x20 = 0.3;
r0 = 0.5;
ctarget = zeros(prod(Ntarget+1), 1);
ctarget(1) = x10^2 + x20^2 - r0^2;
ctarget(2) = -2*x10;
ctarget(6) = -2*x20;
ctarget(3) = 1.0;
ctarget(11) = 1.0;
ctarget(4) = -1.1637;
% ctarget(5) = randn(1);

Btarget = Polynomial.bernstein_transform_matrix(Ntarget);
Etarget = Polynomial.expand_matrix(Ntarget, I);
htarget = prod((1./(Itarget+1)).*( ones(size(Itarget)).^(Itarget+1) - (-ones(size(Itarget))).^(Itarget+1) ),2);

I_q = I( sum(I,2) <= 2, : );
E_q = Polynomial.expand_matrix(N, I_q);
N_lambda = N;
I_lambda = Polynomial.integer_grid(N_lambda);
B_lambda = Polynomial.bernstein_transform_matrix(N_lambda);
h_lambda = prod((1./(I_lambda+1)).*( ones(size(I_lambda)).^(I_lambda+1) - (-ones(size(I_lambda))).^(I_lambda+1) ),2);
P_lambda_x = Polynomial.product_matrix(N_lambda, N);

N_total = N_lambda + N;
E_total = Polynomial.expand_matrix(N_total, I);
E_total_q = Polynomial.expand_matrix(N_total, I_q);
B_total = Polynomial.bernstein_transform_matrix(N_total);

T(length(x1gr)-1, length(x2gr)-1) = struct('gr', [], 'tf', [], 'tf_total', [], 'tf_lambda', [], 'c', [], 'lambda', []);
for i = 1:length(x1gr)-1
    for j = 1:length(x2gr)-1
        
        a_ = 0.5*[x1gr(i+1)-x1gr(i); x2gr(j+1)-x2gr(j)];
        b_ = 0.5*[x1gr(i+1)+x1gr(i); x2gr(j+1)+x2gr(j)];
        
%         T(i,j).tf = Polynomial.domain_transform_matrix(N, a_, b_);
        T(i,j).tf = Polynomial.domain_transform_matrix(Ntarget, a_, b_);
        T(i,j).tf_total = Polynomial.domain_transform_matrix(N_total, a_, b_);
        T(i,j).tf_lambda = Polynomial.domain_transform_matrix(N_lambda, a_, b_);
        
        T(i,j).gr = createGrid( [x1gr(i),x2gr(j)], [x1gr(i+1),x2gr(j+1)], [201, 201]);
    end
end

lambda = zeros( size(I_lambda,1), numel(T) );
%% optimization of V_hat
cvx_begin
    variable c0( size(I,1), 1 )
    variable vi( size(x1Nu{1},2), length(x1Nu) )
    variable vj( size(x1Nu{1},2), length(x2Nu) )
    
    cost = 0.0;
    
    ctmp_i_ = c0;
    for i = 1:size(T,1)
        ctmp_ = ctmp_i_;
        for j = 1:size(T,2)
            Btarget\T(i,j).tf * (ctarget - Etarget*ctmp_) >= 0
            cost = cost - (htarget'*(T(i,j).tf*Etarget*ctmp_));
            if j < size(T,2)
                ctmp_ = ctmp_ + x2Nu{j} * vj(:,j);
            end
        end
        if i < size(T,1)
            ctmp_i_ = ctmp_i_ + x1Nu{i} * vi(:,i);
        end
    end
    
    minimize cost
    
cvx_end

costHist = [];
dxHist = [];
for iter = 1:100
%% optimization of v_tilde
if iter == 1
    cvx_begin
    variable c_tilde( size(I_q,1), 1 )
    
    cost = 0.0;
    cnt = 1;
    ctmp_i_ = c0;
    for i = 1:size(T,1)
        ctmp_ = ctmp_i_;
        for j = 1:size(T,2)
            
            (B_total \ T(i,j).tf_total)*( P_lambda_x*kron(lambda(:,cnt), E_q*c_tilde) + E_total*ctmp_ - E_total_q*c_tilde ) >= 0
            cost = cost - (htarget'*(T(i,j).tf*Etarget*E_q*c_tilde));
            
            if j < size(T,2)
                ctmp_ = ctmp_ + x2Nu{j} * vj(:,j);
            end
            cnt = cnt + 1;
        end
        if i < size(T,1)
            ctmp_i_ = ctmp_i_ + x1Nu{i} * vi(:,i);
        end
    end
    
    minimize cost
    cvx_end
else
    cvx_begin
%     variable c_tilde( size(I_q,1), 1 )
    variable dc_tilde( size(I_q,1), 1 )
    
    c_tilde = c_tilde_prev + dc_tilde;
    
    cnt = 1;
    ctmp_i_ = c0;
    for i = 1:size(T,1)
        ctmp_ = ctmp_i_;
        for j = 1:size(T,2)
            
            (B_total \ T(i,j).tf_total)*( P_lambda_x*kron(lambda(:,cnt), E_q*c_tilde) + E_total*ctmp_ - E_total_q*c_tilde ) >= 0
            
            if j < size(T,2)
                ctmp_ = ctmp_ + x2Nu{j} * vj(:,j);
            end
            cnt = cnt + 1;
        end
        if i < size(T,1)
            ctmp_i_ = ctmp_i_ + x1Nu{i} * vi(:,i);
        end
    end
    
%     norm(dc_tilde) <= 0.5
    
%     minimize dVdc_*c_tilde
    cost = dc_tilde'*(dVdc_'*dVdc_)*dc_tilde + dVdc_*dc_tilde;
    minimize cost
    
    cvx_end
    
%     c_diff = c_tilde - c_tilde_prev;
    c_diff = dc_tilde;
    dxHist = [dxHist, norm(c_diff)];
end


c_tilde_prev = c_tilde;

% compute the next gradient
H_ = [c_tilde(3), 0.5*c_tilde(5); 0.5*c_tilde(5), c_tilde(6)];
b_ = 0.5*c_tilde([2,4]);
r_ = c_tilde(1);

Q_ = (b_'*inv(H_)*b_ - r_)*inv(H_);
invQ_ = inv(Q_);
q_ = -(H_\b_);

S_ = inv(H_);
Svec_ = pinv(Du)*S_(:);

dVdQ_ = invQ_(:)'*Du;
dQdxi2_ = [-Svec_, 2*Svec_*(b_'*S_), Svec_*kron(b_',b_')*Du + (b_'*S_*b_ - r_)*eye(n*(n+1)/2)];
dxi2dxi_ = blkdiag( eye(n+1), -pinv(Du)*kron(S_,S_)*Du );
dVdc_ = dVdQ_ * dQdxi2_ * dxi2dxi_;


costHist = [costHist, log(det(Q_))];


%%
cvx_begin
%     variable lambda( size(I_lambda,1), numel(T) )
    variable lambda0( size(I,1), 1 )
    variable lambda_i( size(x1Nu{1},2), length(x1Nu) )
    variable lambda_j( size(x1Nu{1},2), length(x2Nu) )
    
%     cost = 0.0;
    cnt = 1;
    ctmp_i_ = c0;
    lambdatmp_i_ = lambda0;
    
    for i = 1:size(T,1)
        ctmp_ = ctmp_i_;
        lambdatmp_ = lambdatmp_i_;
        
        for j = 1:size(T,2)
            
%             (B_total \ T(i,j).tf_total)*( P_lambda_x*kron(lambda(:,cnt), E_q*c_tilde) + E_total*ctmp_ - E_total_q*c_tilde ) >= 0
%             (B_lambda \ T(i,j).tf_lambda)* lambda(:,cnt) >= 0
            (B_total \ T(i,j).tf_total)*( P_lambda_x*kron(lambdatmp_, E_q*c_tilde) + E_total*ctmp_ - E_total_q*c_tilde ) >= 0
            (B_lambda \ T(i,j).tf_lambda)* lambdatmp_ >= 0
%             cost = cost + (h_lambda'*(T(i,j).tf_lambda*lambda(:,cnt)));
            
            if j < size(T,2)
                ctmp_ = ctmp_ + x2Nu{j} * vj(:,j);
                lambdatmp_ = lambdatmp_ + x2Nu{j} * lambda_j(:,j);
            end
            cnt = cnt + 1;
        end
        if i < size(T,1)
            ctmp_i_ = ctmp_i_ + x1Nu{i} * vi(:,i);
            lambdatmp_i_ = lambdatmp_i_ + x1Nu{i} * lambda_i(:,i);
        end
    end

%     minimize cost
cvx_end

%%
% c0 = randn(size(I,1), 1);
% vi1 = randn( size(x1Nu{1},2), 1 );
% vi2 = randn( size(x1Nu{2},2), 1 );
% vj1 = randn( size(x2Nu{1},2), 1 );
% vj2 = randn( size(x2Nu{2},2), 1 );

cnt = 1;
ctmp_i_ = c0;
lambdatmp_i_ = lambda0;
for i = 1:size(T,1)
    ctmp_ = ctmp_i_;
    lambdatmp_ = lambdatmp_i_;
    for j = 1:size(T,2)
        T(i,j).c = ctmp_;
        T(i,j).lambda = lambdatmp_;
        lambda(:,cnt) = lambdatmp_;
%         T(i,j).lambda = lambda(:,cnt);
        if j < size(T,2)
            ctmp_ = ctmp_ + x2Nu{j} * vj(:,j);
            lambdatmp_ = lambdatmp_ + x2Nu{j} * lambda_j(:,j);
        end
        cnt = cnt + 1;
    end
    if i < size(T,1)
        ctmp_i_ = ctmp_i_ + x1Nu{i} * vi(:,i);
        lambdatmp_i_ = lambdatmp_i_ + x1Nu{i} * lambda_i(:,i);
    end
end

gr = createGrid([-1,-1],[1,1],[51,51]);
ptarget = Polynomial( ctarget, Itarget );

p_q = Polynomial( c_tilde, I_q );

figure(321)
cla; hold on; grid on;
for i = 1:size(T,1)
    for j = 1:size(T,2)
        p_ = Polynomial( T(i,j).c, I );
        surf( T(i,j).gr.xs{1}, T(i,j).gr.xs{2}, p_.eval( T(i,j).gr ), 'linestyle', 'none', 'facealpha', 0.5, 'facecolor', rand(3,1))
    end
end
surf(gr.xs{1}, gr.xs{2}, ptarget.eval( gr ), 'linestyle', '-', 'facealpha', 0.2, 'facecolor', 'none')
surf(gr.xs{1}, gr.xs{2}, p_q.eval( gr ), 'linestyle', '-', 'facealpha', 0.2, 'facecolor', 'none', 'edgecolor', 'b')
xlabel('$x_1$')
ylabel('$x_2$')
title('$\hat{V}(x)$ and $\tilde{V}(x)$')
view([0,0])
% zlim([-0.01, 0.06])

figure(325)
cla; hold on; grid on;
for i = 1:size(T,1)
    for j = 1:size(T,2)
%         p_ = Polynomial( T(i,j).c, I );
        p_ = Polynomial( T(i,j).lambda, I_lambda );
        surf( T(i,j).gr.xs{1}, T(i,j).gr.xs{2}, p_.eval( T(i,j).gr ), 'linestyle', 'none', 'facealpha', 0.5, 'facecolor', rand(3,1))
    end
end
% surf(gr.xs{1}, gr.xs{2}, ptarget.eval( gr ), 'linestyle', '-', 'facealpha', 0.2, 'facecolor', 'none')
% surf(gr.xs{1}, gr.xs{2}, p_q.eval( gr ), 'linestyle', '-', 'facealpha', 0.2, 'facecolor', 'none', 'edgecolor', 'b')
xlabel('$\lambda_1$')
ylabel('$\lambda_2$')
title('$\Lambda(x)$')
view([0,0])

% figure(322)
% cla; hold on; grid on;
% surf(gr.xs{1}, gr.xs{2}, p_q.eval( gr ), 'linestyle', '-', 'facealpha', 0.2, 'facecolor', 'none', 'edgecolor', 'b')
% xlabel('$x_1$')
% ylabel('$x_2$')
% view([0,0])

figure(323)
cla; hold on; grid on; axis equal;
tmp = Utils.get_level_set( gr, ptarget.eval( gr ), 0.0 );
plot(tmp(1,:), tmp(2,:), 'k-', 'linewidth', 2)
% tmp = Utils.get_level_set( gr, p_q.eval( gr ), 0.0 );
tmp = q_ + Q_^(1/2) * Math.Sphere(1,200).x;
plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2)
xlabel('$x_1$')
ylabel('$x_2$')
xlim([-1,1])
ylim([-1,1])
drawnow

figure(328)
cla; hold on; grid on; axis tight;
plot(1:iter, costHist, 'bx-', 'linewidth', 2)
xlabel('Iteration')
ylabel('Cost')
drawnow

% if iter > 1
%     figure(329)
%     cla; hold on; grid on; axis tight;
%     plot(1:(iter-1), dxHist, 'bx-', 'linewidth', 2)
%     xlabel('Iteration')
%     ylabel('Rate')
%     drawnow
% end
end
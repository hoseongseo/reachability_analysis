clear; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

N = 4;
I = Polynomial.integer_grid(N);
% c = randn( prod(N+1), 1 );
c = [-0.7481; 0.2046; 0.2750; 1.6998; 1.6156];
p = Polynomial(c, I );

B = inv(Polynomial.bernstein_transform_matrix(N));
h = prod((1./(I+1)).*( ones(size(I)).^(I+1) - (-ones(size(I))).^(I+1) ),2);


N_hat = 1;
I_hat = Polynomial.integer_grid(N_hat);
h_hat = prod((1./(I_hat+1)).*( ones(size(I_hat)).^(I_hat+1) - (-ones(size(I_hat))).^(I_hat+1) ),2);
E_hat = Polynomial.expand_matrix(N, I_hat);

%%
% without subdivision
cvx_begin
    variable c_hat1( prod(N_hat+1), 1 )
    
    minimize -(h_hat'*c_hat1)
%     minimize max( B*(c-E_hat*c_hat1) )
    subject to
    B * (c - E_hat * c_hat1) >= 0

cvx_end
p_hat1 = Polynomial(c_hat1, I_hat);
disp(['cost = ', num2str(-h_hat'*c_hat1)])

d = B*c;
C = B*E_hat;
% [c2,c3] = meshgrid(linspace(-3,3,101), linspace(-3,3,101));
c2 = [-3;3];

figure(1)
cla; hold on; grid on; axis tight; axis equal;
xlabel('$c_1$')
ylabel('$c_2$')
pts = [-3,-3; -3,3; 3,-3; 3,3];
for i = 1:size(B,1)
    c1 = (1/C(i,1))*(-C(i,2)*c2 + d(i));
%     plot(c1, c2, 'kx-', 'linewidth', 1)
    
    vals = C(i,1)*pts(:,1) + C(i,2)*pts(:,2);
    idx = vals > d(i);
    
    pts_valid = [c1, c2; pts(idx,:)];
    K = convhull(pts_valid(:,1), pts_valid(:,2));
    patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.3, 'facecolor', 'k', 'linestyle', 'none')
    pause
end
plot(c_hat1(1), c_hat1(2), 'r*')
xlim([-3,0])
ylim([-2,1])


%%
% with subdivision (sparse)
T1 = Polynomial.domain_transform_matrix( N, 0.5, -0.5);
T2 = Polynomial.domain_transform_matrix( N, 0.5, 0.5);
cvx_begin
    variable c_hat2( prod(N_hat+1), 1 )
    
    minimize -(h_hat'*c_hat2)
%     minimize max( B*(c-E_hat*c_hat2) )
    
    subject to
    B * T1 * (c - E_hat * c_hat2) >= 0
    B * T2 * (c - E_hat * c_hat2) >= 0
cvx_end
p_hat2 = Polynomial(c_hat2, I_hat);


c2 = [-3;3];
pts = [-3,-3; -3,3; 3,-3; 3,3];


figure(2)
cla; hold on; grid on; axis tight; axis equal;
xlabel('$c_1$')
ylabel('$c_2$')
% xlim([-3,3])
% ylim([-3,3])
xlim([-3,0])
ylim([-2,1])

d = B*T1*c;
C = B*T1*E_hat;

for i = 1:size(B,1)
    c1 = (1/C(i,1))*(-C(i,2)*c2 + d(i));
%     plot(c1, c2, 'kx-', 'linewidth', 1)
    
    vals = C(i,1)*pts(:,1) + C(i,2)*pts(:,2);
    idx = vals > d(i);
    
    pts_valid = [c1, c2; pts(idx,:)];
    K = convhull(pts_valid(:,1), pts_valid(:,2));
    patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.15, 'facecolor', 'k', 'linestyle', 'none')
end

d = B*T2*c;
C = B*T2*E_hat;
for i = 1:size(B,1)
    c1 = (1/C(i,1))*(-C(i,2)*c2 + d(i));
%     plot(c1, c2, 'kx-', 'linewidth', 1)
    
    vals = C(i,1)*pts(:,1) + C(i,2)*pts(:,2);
    idx = vals > d(i);
    
    pts_valid = [c1, c2; pts(idx,:)];
    K = convhull(pts_valid(:,1), pts_valid(:,2));
    patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.15, 'facecolor', 'k', 'linestyle', 'none')
end
plot(c_hat1(1), c_hat1(2), 'r*')
plot(c_hat2(1), c_hat2(2), 'b*')

%%
% with subdivision (fine)
T1 = Polynomial.domain_transform_matrix( N, 0.25, 0.25*-3);
T2 = Polynomial.domain_transform_matrix( N, 0.25, 0.25*-1);
T3 = Polynomial.domain_transform_matrix( N, 0.25, 0.25*1);
T4 = Polynomial.domain_transform_matrix( N, 0.25, 0.25*3);
cvx_begin
    variable c_hat3( prod(N_hat+1), 1 )
    
    minimize -(h_hat'*c_hat3)
%     minimize max( B*(c-E_hat*c_hat3) )

    subject to
    B * T1 * (c - E_hat * c_hat3) >= 0
    B * T2 * (c - E_hat * c_hat3) >= 0
    B * T3 * (c - E_hat * c_hat3) >= 0
    B * T4 * (c - E_hat * c_hat3) >= 0
cvx_end
p_hat3 = Polynomial(c_hat3, I_hat);


figure(3)
cla; hold on; grid on; axis tight; axis equal;
xlabel('$c_1$')
ylabel('$c_2$')
% xlim([-3,3])
% ylim([-3,3])
% xlim([-3,0])
% ylim([-2,1])
xlim([-0.85,-0.65])
ylim([0.2, 0.4])

d = B*T1*c;
C = B*T1*E_hat;
for i = 1:size(B,1)
    c1 = (1/C(i,1))*(-C(i,2)*c2 + d(i));
%     plot(c1, c2, 'kx-', 'linewidth', 1)
    
    vals = C(i,1)*pts(:,1) + C(i,2)*pts(:,2);
    idx = vals > d(i);
    
    pts_valid = [c1, c2; pts(idx,:)];
    K = convhull(pts_valid(:,1), pts_valid(:,2));
    patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.1, 'facecolor', 'k', 'linestyle', 'none')
end

d = B*T2*c;
C = B*T2*E_hat;
for i = 1:size(B,1)
    c1 = (1/C(i,1))*(-C(i,2)*c2 + d(i));
%     plot(c1, c2, 'kx-', 'linewidth', 1)
    
    vals = C(i,1)*pts(:,1) + C(i,2)*pts(:,2);
    idx = vals > d(i);
    
    pts_valid = [c1, c2; pts(idx,:)];
    K = convhull(pts_valid(:,1), pts_valid(:,2));
    patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.1, 'facecolor', 'k', 'linestyle', 'none')
end

d = B*T3*c;
C = B*T3*E_hat;
for i = 1:size(B,1)
    c1 = (1/C(i,1))*(-C(i,2)*c2 + d(i));
%     plot(c1, c2, 'kx-', 'linewidth', 1)
    
    vals = C(i,1)*pts(:,1) + C(i,2)*pts(:,2);
    idx = vals > d(i);
    
    pts_valid = [c1, c2; pts(idx,:)];
    K = convhull(pts_valid(:,1), pts_valid(:,2));
    patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.1, 'facecolor', 'k', 'linestyle', 'none')
end

d = B*T4*c;
C = B*T4*E_hat;
for i = 1:size(B,1)
    c1 = (1/C(i,1))*(-C(i,2)*c2 + d(i));
%     plot(c1, c2, 'kx-', 'linewidth', 1)
    
    vals = C(i,1)*pts(:,1) + C(i,2)*pts(:,2);
    idx = vals > d(i);
    
    pts_valid = [c1, c2; pts(idx,:)];
    K = convhull(pts_valid(:,1), pts_valid(:,2));
    patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.1, 'facecolor', 'k', 'linestyle', 'none')
end

plot(c_hat1(1), c_hat1(2), 'r*')
plot(c_hat2(1), c_hat2(2), 'b*')
plot(c_hat3(1), c_hat3(2), '*', 'color', [39,176,39]/255)

%%
% another approach
a = sym('a', 'real'); % domain transform parameters
b = sym('b', 'real');
T = Polynomial.domain_transform_matrix( N, a, b );

rhs = B*T*c;
c_hat = sym('c_hat', [prod(N_hat+1),1], 'real');
lhs = B*T*E_hat*c_hat;

% lhs_b = subs(lhs, a, eps);
% rhs_b = subs(rhs, a, eps);
lhs_b = subs(lhs, a, 0.0);
rhs_b = subs(rhs, a, 0.0);
geq = rhs_b - lhs_b;

% %%
% geq1 = geq(1);
% geq1_poly = Polynomial.fromSym( geq1, b );
% N_poly = max(geq1_poly.order,[],1);
% geq1_b_coeff = Polynomial.bernstein_transform_matrix(N_poly) * geq1_poly.coeff; % this must be greater than 0
% geq1_A = jacobian(geq1_b_coeff, c_hat);
% geq1_B = geq1_b_coeff - geq1_A * c_hat;
% geq1_A = double(geq1_A);
% geq1_B = double(geq1_B);

figure(55)
cla; hold on; grid on; axis equal; axis tight;
xlabel('$c_1$')
ylabel('$c_2$')
% xlim([-3,3])
% ylim([-3,3])
% xlim([-3,0])
% ylim([-2,1])
xlim([-0.85,-0.65])
ylim([0.2, 0.4])
for i = 1%:length(geq)
    geq_ = geq(i);
    for b_test = linspace(-1,1,101)
        geq1_ = subs(geq_, b, b_test);
        geq1_A_ = jacobian(geq1_, c_hat);
        geq1_B_ = geq1_ - geq1_A_*c_hat;
        geq1_A_ = double(geq1_A_);
        geq1_B_ = double(geq1_B_);
        
        c1 = (1/geq1_A_(1))*(-geq1_A_(2)*c2 - geq1_B_);
%         plot(c1, c2, 'kx-', 'linewidth', 1)
        
        vals = geq1_A_(1)*pts(:,1) + geq1_A_(2)*pts(:,2);
        idx = vals < -geq1_B_;
        
        pts_valid = [c1, c2; pts(idx,:)];
        K = convhull(pts_valid(:,1), pts_valid(:,2));
        
        alpha_val = 0.05*abs(b_test + 0.3);
        if i == 1
            patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', alpha_val, 'facecolor', 'k', 'linestyle', 'none', 'edgecolor', 'k')
        elseif i == 2
            patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.02, 'facecolor', 'b', 'linestyle', '-', 'edgecolor', 'b')
        elseif i == 3
            patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.02, 'facecolor', 'r', 'linestyle', '-', 'edgecolor', 'r')
        else
            patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.02, 'facecolor', 'g', 'linestyle', '-', 'edgecolor', 'g')
        end
        
%         if (b_test == -1) || (b_test == 1)
%             patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.2, 'facecolor', 'k', 'linestyle', 'none', 'edgecolor', 'k')
%         end 
%         pause
        drawnow
    end
end
plot(c_hat1(1), c_hat1(2), 'r*')
plot(c_hat2(1), c_hat2(2), 'b*')
plot(c_hat3(1), c_hat3(2), '*', 'color', [39,176,39]/255)

%%
% f = @(b) c(1) + b*c(2) + b.^2*c(3) + b.^3*c(4) + b.^4*c(5);
% for b_test = linspace(-1,1,101)
%     c1 = f(b_test) - b_test*c2;
%     plot( c1, c2, 'm--', 'linewidth', 1)
% end
% for b_test = [-0.5, 0.5]
%     c1 = f(b_test) - b_test*c2;
%     plot( c1, c2, 'g--', 'linewidth', 2)
% end

% b_gr = linspace(-1,1,5)';
% A_gr = [ones(size(b_gr)), b_gr];
% B_gr = f(b_gr);
% cvx_begin
%     variable c_gr( prod(N_hat+1), 1 )
%     minimize -(h_hat'*c_gr)
%     subject to
%     A_gr * c_gr <= B_gr
% cvx_end
% p_gr = Polynomial(c_gr, I_hat);
% 
% plot(c_gr(1), c_gr(2), 'y*')
%%

% %%
% figure(44)
% cla; hold on; grid on; axis equal; axis tight;
% xlabel('$c_1$')
% ylabel('$c_2$')
% xlim([-3,3])
% ylim([-3,3])
% for i = 1:size(B,1)
%     c1 = (1/geq1_A(i,1))*(-geq1_A(i,2)*c2 - geq1_B(i));
%     plot(c1, c2, 'kx--', 'linewidth', 2)
%     
%     vals = geq1_A(i,1)*pts(:,1) + geq1_A(i,2)*pts(:,2);
%     idx = vals < -geq1_B(i);
%     
%     pts_valid = [c1, c2; pts(idx,:)];
%     K = convhull(pts_valid(:,1), pts_valid(:,2));
%     patch('xdata', pts_valid(K,1), 'ydata', pts_valid(K,2), 'facealpha', 0.2, 'facecolor', 'k')
% end
% 
% plot(c_hat1(1), c_hat1(2), 'r*')
% plot(c_hat2(1), c_hat2(2), 'b*')
% plot(c_hat3(1), c_hat3(2), 'go')



%%
xgr = createGrid(-1, 1, 101);

figure(11)
cla; hold on; grid on; axis tight;


plot(xgr.xs{1}, p_hat1.eval(xgr), 'r', 'linewidth', 2)
plot(xgr.xs{1}, p_hat2.eval(xgr), 'b-', 'linewidth', 2)
plot(xgr.xs{1}, p_hat3.eval(xgr), 'color', [39,176,39]/255, 'linewidth', 3)
plot(xgr.xs{1}, p.eval(xgr), 'k', 'linewidth', 2)

xlabel('$x$')
ylabel('$V$')
xlim([-1,1])
ylim([-3.8589, 3.0469])
% ylim([0,1])

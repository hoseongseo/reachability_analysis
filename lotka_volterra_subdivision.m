clear; clc; 

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

% given initial set
Q0 = diag([0.05; 0.05])^2;

% disturbance bound
wMax = 0.1;

% sphere
S = Utils.Sphere(1,200);

%%% system (polynomial dynamics)
t = linspace(0,1,51);
q0 = [1.2; 1.1];
sys0 = Dynamics.LotkaVolterraNominal([1.2; 1.1], t); % nominal dynamics with the initial condition
% sysN = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin
sys = Dynamics.LotkaVolterra(sys0, wMax); % system shifted to the origin

% sys = Dynamics.LotkaVolterra2( wMax );

%% subdivision
N = [2,2];

% declare state and disturbance
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

N_a = [N, zeros(1,sys.Nw)]; % augmentation
fp = sys.f(xp, zeros(sys.Nu,1), wp, t(1));
fp_order_all = [];
for i = 1:sys.Nx
    fp_order_all = [fp_order_all; fp(i).order];
end
N_f = max(reshape(fp_order_all, [numel(fp_order_all)/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);
N_max = N_a + N_f; % total order of polynmial approximation 


E_N = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a));
E_N_N_f = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a + N_f));

B_N = Polynomial.bernstein_transform_matrix(N);
B_N_max = Polynomial.bernstein_transform_matrix(N_max);

P_N_N_f = Polynomial.product_matrix(N_a, N_f);

D = Polynomial.derivative_matrix(N);

M = zeros( prod(N_a+N_f+1), prod(N+1), length(t)-1 );
for k = 1:length(t)-1
    mat = zeros(prod(N_a+N_f+1), prod(N+1));
%     if isfield(args, 'uN')
%         f_ = sys.f(xp, args.uN(:,k), wp, t(k));
%     else
        f_ = sys.f(xp, zeros(sys.Nu,1), wp, t(k));
%     end
    for j = 1:sys.Nx
        tmp = P_N_N_f * kron( eye(prod(N+1)), Polynomial.expand_matrix(N_f, f_(j).order)*f_(j).coeff ) * D(:,:,j);
        mat = mat + tmp;
    end
    M(:,:,k) = mat;
end

poly_degree = Polynomial.integer_grid(N);
quad_degree = poly_degree((sum(poly_degree, 2) == 2) | (sum(poly_degree, 2) == 0), :);
% quad_degree = poly_degree( (sum(poly_degree, 2) <= 2), :);
h_poly = prod((1./(poly_degree+1)).*( ones(size(poly_degree)).^(poly_degree+1) - (-ones(size(poly_degree))).^(poly_degree+1) ),2);

E_quad_poly = Polynomial.expand_matrix(N, quad_degree);

[~, full_order_quad] = Polynomial.affine_transform_matrix2(quad_degree, Q0^(1/2), zeros(sys.Nx,1));
h_quad = prod((1./(full_order_quad+1)).*( ones(size(full_order_quad)).^(full_order_quad+1) - (-ones(size(full_order_quad))).^(full_order_quad+1) ),2);

Q_ltv = funnel_ltv(sys, t, Q0);
% region = struct('q', sys0.xN, 'Q', Q_ltv * 1.2);
region = struct('q', zeros(size(sys0.xN)), 'Q', Q_ltv * 1.2);

V0 = xp'*inv(Q0)*xp - 1;
% V0 = (xp-q0)'*inv(Q0)*(xp-q0) - 1;
c0 = Polynomial.expand_matrix(N, V0.order(:,1:sys.Nx))*V0.coeff;
c_quad0 = c0((sum(poly_degree, 2) == 2) | (sum(poly_degree, 2) == 0));
% c_quad0 = c0(sum(poly_degree, 2) <= 2);

% %%
% Vtmp = Polynomial( c0, poly_degree );
% gr = createGrid(q0 + [-0.2;-0.2], q0 + [0.2;0.2], [201;201]);

% figure(11)
% cla; hold on; grid on; axis tight;
% surf(gr.xs{1}, gr.xs{2}, Vtmp.eval( gr ) )

%% subdivision parameter
n_div_max = [2; 2; 1];
argin = cell(1,length(n_div_max));
for j = 1:length(n_div_max)
    argin{j} = 1:n_div_max(j);
end
argout = cell(1,length(n_div_max));
[argout{:}] = ndgrid(argin{:});
n_div_max_gr = zeros(prod(n_div_max),length(n_div_max));
for j = 1:length(n_div_max)
    n_div_max_gr(:,j) = argout{j}(:);
end

n_div = [2; 2];
argin = cell(1,length(n_div));
for j = 1:length(n_div)
    argin{j} = 1:n_div(j);
end
argout = cell(1,length(n_div));
[argout{:}] = ndgrid(argin{:});
n_div_gr = zeros(prod(n_div),length(n_div));
for j = 1:length(n_div)
    n_div_gr(:,j) = argout{j}(:);
end

region_ = region(end);

T_N_max = zeros( prod(N_max+1), prod(N_max+1), length(t), prod(n_div_max));
T_N = zeros( prod(N+1), prod(N+1), length(t), prod(n_div) );
for k = 1:length(t)
    sqrtQ_ = region_.Q(:,:,k)^(1/2);
    a_ = sqrt(sum(sqrtQ_.^2,2));
    b_ = region_.q(:,k);
    for i = 1:prod(n_div_max)
        n_div_max_gr_ = n_div_max_gr(i,:)';
        T_N_max(:,:,k,i) = Polynomial.domain_transform_matrix( N_max,...
            (1./n_div_max).*[a_; ones(sys.Nw,1)],...
            [b_; zeros(sys.Nw,1)] + (1./n_div_max).*[a_; ones(sys.Nw,1)].*( 2*(n_div_max_gr_-1) - (n_div_max-1) ) );
    end
    
    for i = 1:prod(n_div)
        n_div_gr_ = n_div_gr(i,:)';
        T_N(:,:,k,i) = Polynomial.domain_transform_matrix( N,...
            (1./n_div).*a_,...
            (1./n_div).*a_.*( 2*(n_div_gr_-1) - (n_div-1) ) );
    end
end



cvx_begin
    variable c( prod(N+1), length(t) )
    variable c_quad( (sys.Nx)*(sys.Nx+1)/2+1, length(t) )
    
%     variable c( length(c_quad0), length(t) )
    cost = 0.0;
    
    % IC
    c(:,1) == c0
    c_quad(:,1) == c_quad0
%     c(:,1) == c_quad0
    
    for k = 1:length(t)
        sqrtQ_ = region_.Q(:,:,k)^(1/2);
        a_ = sqrt(sum(sqrtQ_.^2,2));
        b_ = region_.q(:,k);
        
        % CONSERVATIVENESS
        if k < length(t)
            for i = 1:prod(n_div_max)
%                 n_div_max_gr_ = n_div_max_gr(i,:)';
%                 T_N_max_ = Polynomial.domain_transform_matrix( N_max,...
%                     (1./n_div_max).*[a_; ones(sys.Nw,1)],...
%                     [b_; zeros(sys.Nw,1)] + (1./n_div_max).*[a_; ones(sys.Nw,1)].*( 2*(n_div_max_gr_-1) - (n_div_max-1) ) );
%                 B_N_max\T_N_max_*( E_N*(c(:,k+1)-c(:,k)) + (t(k+1)-t(k))*(E_N_N_f*M(:,:,k)*c(:,k)) ) <= 0
                B_N_max\T_N_max(:,:,k,i)*( E_N*(c(:,k+1)-c(:,k)) + (t(k+1)-t(k))*(E_N_N_f*M(:,:,k)*c(:,k)) ) <= 0
            end
        end
        if k > 1
            for i = 1:prod(n_div)
%                 n_div_gr_ = n_div_gr(i,:)';
%                 T_N_ = Polynomial.domain_transform_matrix( N,...
%                     (1./n_div).*a_,...
%                     (1./n_div).*a_.*( 2*(n_div_gr_-1) - (n_div-1) ) );
%                 B_N\T_N_*( E_quad_poly*c_quad(:,k) - c(:,k) ) <= 0
                B_N\T_N(:,:,k,i)*( E_quad_poly*c_quad(:,k) - c(:,k) ) <= 0
            end
        end
        
        % COST 
        if k > 1
            T_quad_ = Polynomial.affine_transform_matrix2( quad_degree, sqrtQ_, b_ );
            cost = cost - h_quad'*(T_quad_*c_quad(:,k));
        end
    end
    
    minimize cost
cvx_end

%% HJB
nGrid = [201; 201];
minGrid = [-0.2; -0.2];
maxGrid = [0.2; 0.2];
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
% figure(221)
figure;
cla; hold on; grid on; axis tight; axis equal;
for k = round(linspace(1,length(t),5))
%     tmp = sys0.xN(:,k) + Q_ltv(:,:,k)^(1/2) * Math.Sphere(1,200).x;
%     plot(tmp(1,:), tmp(2,:), 'g-', 'linewidth', 2)
    
    tmp = sys0.xN(:,k) + X{k};
    plot(tmp(1,:), tmp(2,:), 'k-', 'linewidth', 2)
    
    c_quad_ = c_quad(:,k);
    Q_ = (-c_quad_(1))*inv([c_quad_(2), 0.5*c_quad_(3); 0.5*c_quad_(3), c_quad_(4)]);
    tmp = sys0.xN(:,k) + Q_^(1/2) * Math.Sphere(1,200).x;
    plot(tmp(1,:), tmp(2,:), 'b-', 'linewidth', 2)
    
    
    Vp = Polynomial(c(:,k), poly_degree);
    gr_ = createGrid(minGrid, maxGrid, nGrid);
    tmp = sys0.xN(:,k) + Utils.get_level_set( gr_, Vp.eval(gr_), 0.0 );
    if k == 1
        plot(tmp(1,:), tmp(2,:), 'k', 'linewidth', 2)
    else
        plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2)
    end
    
end

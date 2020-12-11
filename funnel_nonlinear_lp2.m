function [region, cost_hist, rate_hist, c, s] = funnel_nonlinear_lp2(sys, tk, Q0, N, args)
% sys: dynamics
% tk: time horizon of interest (discretized)
% Q0: shape matrix for the initial set
% N: target order of the approximated value function
% args: extra arguments detailed as
if nargin < 5
    args.max_iter = 1;
    args.ftol = 1e-5;
    args.plot_cost = false;
end

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

% domain of interest (initialized by the funnel of LTV system)
Q_ltv = funnel_ltv(sys, tk, Q0);
region = struct('q', zeros(sys.Nx,length(tk)), 'Q', Q_ltv);

% degrees of the approximated value functions
N_aug = [N, zeros(1,sys.Nw)];
fp = sys.f(xp, zeros(sys.Nu,1), wp, tk(1));
N_dyn = max(reshape([fp.order], [numel([fp.order])/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);
% N_dyn = [2,2,2,2];

N_all = N_aug + N_dyn; % total order of approximation

E = Polynomial.expand_matrix(N_all, Polynomial.integer_grid(N_aug));
B_N_all = Polynomial.bernstein_transform_matrix( N_all );
P = Polynomial.product_matrix(N_aug, N_dyn);
B_N_dyn = Polynomial.bernstein_transform_matrix( N_dyn );
D = Polynomial.derivative_matrix(N_aug);
B_N = Polynomial.bernstein_transform_matrix( N );

poly_degree = Polynomial.integer_grid(N);
% quad_degree = poly_degree(sum(poly_degree, 2) <= 2, :);
quad_degree = poly_degree((sum(poly_degree, 2) == 2) | (sum(poly_degree, 2) == 0), :);
% quad_idx_1st_order = sum(quad_degree,2) == 1;
quad_idx_2nd_order = sum(quad_degree,2) == 2;
quad_idx_skew = max(quad_degree,[],2) == 1;
E_N = Polynomial.expand_matrix( N, quad_degree );

% initial coefficients
V0 = xp'*inv(Q0)*xp - 1;
c0 = Polynomial.expand_matrix(N, V0.order(:,1:sys.Nx))*V0.coeff;
% c_tilde0 = c0(sum(poly_degree, 2) <= 2);
c_tilde0 = c0((sum(poly_degree, 2) == 2) | (sum(poly_degree, 2) == 0));

%%%%%%%%%
% % order of dynamics
% fp = sys.f(xp, zeros(sys.Nu,1), wp, tk(1));
% N_dyn = max(reshape([fp.order], [numel([fp.order])/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);
% Nx = N_dyn(1:sys.Nx); % order in state
% Nw = N_dyn(sys.Nx + (1:sys.Nw)); % order in disturbance
% 
% % precomputation of matrices
% D_N = Polynomial.derivative_matrix( N );
% invB_Nw = inv(Polynomial.bernstein_transform_matrix( Nw ));
% B_N = Polynomial.bernstein_transform_matrix( N );
% 
% P_N_Nx = Polynomial.product_matrix(N, Nx);
% I_N = eye(prod(N+1));
% C_N_Nx = Polynomial.expand_matrix( N+Nx, Polynomial.integer_grid(N));
% E_N = Polynomial.expand_matrix( N, quad_degree );
% 
% B_NNx = (Polynomial.bernstein_transform_matrix( N+Nx ));
h_N = prod((1./(poly_degree+1)).*( ones(size(poly_degree)).^(poly_degree+1) - (-ones(size(poly_degree))).^(poly_degree+1) ),2);

% Hamiltonian-bounding matrices over time
M = zeros( prod(N_all+1), prod(N+1), length(tk)-1 );
for k = 1:length(tk)-1
    mat = zeros(prod(N_all+1), prod(N+1));
    f_ = sys.f(xp, zeros(sys.Nu,1), wp, tk(k));
    for j = 1:sys.Nx
        mat = mat + P * kron( eye(prod(N_aug+1)), Polynomial.expand_matrix(N_dyn, f_(j).order)*f_(j).coeff ) * D(:,:,j);
    end
    M(:,:,k) = mat;
end


lambda = zeros( prod(N_dyn+1), length(tk)-1 ); % space for lambda
lambda_tilde = zeros( prod(N+1), length(tk)-1 ); 

P_2N = Polynomial.product_matrix(N, N);
B_2N = Polynomial.bernstein_transform_matrix( 2*N );
T_2N = zeros( prod(2*N+1), prod(2*N+1), length(tk) );
E_2N = Polynomial.expand_matrix( 2*N, quad_degree );

E_2N_poly = Polynomial.expand_matrix( 2*N, poly_degree );

% volume of the initial guess
vol = 0.0;
for k = 2:length(tk)
%     vol = vol + log( area(polyshape(  (region.Q(:,:,k)^(1/2)*Utils.Sphere(1,800).x)'  )) );
    vol = vol + log(det( region.Q(:,:,k) ));
%     vol = vol + trace( region.Q(:,:,k) );
end
cost_hist = vol;
rate_hist = [];


T_N_all = zeros( prod(N_all+1), prod(N_all+1), length(tk) );
T_N_dyn = zeros( prod(N_dyn+1), prod(N_dyn+1), length(tk) );
T_N = zeros( prod(N+1), prod(N+1), length(tk) );

[~, full_order] = Polynomial.affine_transform_matrix(N, Q0^(1/2), zeros(sys.Nx,1));
T_N2 = zeros( size(full_order,1), prod(N+1), length(tk) );
h_N2 = prod((1./(full_order+1)).*( ones(size(full_order)).^(full_order+1) - (-ones(size(full_order))).^(full_order+1) ),2);

[~, full_order_tilde] = Polynomial.affine_transform_matrix2(quad_degree, Q0^(1/2), zeros(sys.Nx,1));
T_N_tilde = zeros( size(full_order_tilde,1), size(quad_degree,1), length(tk) );
h_N_tilde = prod((1./(full_order_tilde+1)).*( ones(size(full_order_tilde)).^(full_order_tilde+1) - (-ones(size(full_order_tilde))).^(full_order_tilde+1) ),2);

    region_ = region(end);
    for k = 1:length(tk)
        sqrtQ_ = region_.Q(:,:,k)^(1/2);
        a_ = sqrt(sum(sqrtQ_.^2,2));
        b_ = region_.q(:,k);
    
        T_N_all(:,:,k) = Polynomial.domain_transform_matrix( N_all, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
        T_N_dyn(:,:,k) = Polynomial.domain_transform_matrix( N_dyn, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
        T_N(:,:,k) = Polynomial.domain_transform_matrix( N, a_, b_ );
        T_2N(:,:,k) = Polynomial.domain_transform_matrix( 2*N, a_, b_ );
%         T_N_all(:,:,k) = eye(prod(N_all+1));
%         T_N_dyn(:,:,k) = eye(prod(N_dyn+1));
%         T_N(:,:,k) = eye(prod(N+1));
%         T_2N(:,:,k) = eye(prod(2*N+1));
        T_N2(:,:,k) = Polynomial.affine_transform_matrix( N, sqrtQ_, b_ );
        T_N_tilde(:,:,k) = Polynomial.affine_transform_matrix2( quad_degree, sqrtQ_, b_ );
    end


% MAIN LOOP
num_elem = 21; % number of elements per segment
num_seg = (length(tk)-1)/(num_elem-1);

for iter = 1:args.max_iter
%     %% Domain transform matrices over time
    % current domain of interest
    
%     region_ = region(end);
%     for k = 1:length(tk)
%         sqrtQ_ = region_.Q(:,:,k)^(1/2);
%         a_ = sqrt(sum(sqrtQ_.^2,2));
%         b_ = region_.q(:,k);
%     
%         T_N_all(:,:,k) = Polynomial.domain_transform_matrix( N_all, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
%         T_N_dyn(:,:,k) = Polynomial.domain_transform_matrix( N_dyn, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
%         T_N(:,:,k) = Polynomial.domain_transform_matrix( N, a_, b_ );
%         T_2N(:,:,k) = Polynomial.domain_transform_matrix( 2*N, a_, b_ );
%         %     T_N2(:,:,k) = Polynomial.affine_transform_matrix( N, sqrtQ_, b_ );
%     end
    
    %% Solve incrementally
%     if iter == 1
        c = zeros(prod(N+1), length(tk));
        c_tilde = zeros((sys.Nx)*(sys.Nx+1)/2+1, length(tk));
                
        for i = 1:num_seg
            cvx_begin
                cvx_precision best
%                 cvx_quiet true
            
                variable c_seg( prod(N+1), num_elem )
                variable c_tilde_seg( (sys.Nx)*(sys.Nx+1)/2+1, num_elem )
                
                cost = 0.0;
                for k = 2:num_elem
                    idx_k = (num_elem-1)*(i-1)+k;
%                     cost = cost -h_N'*(T_N(:,:,idx_k)*c_seg(:,k));
%                     cost = cost -h_N2'*(T_N2(:,:,idx_k)*c_seg(:,k));
%                     cost = cost -h_N'*(T_N(:,:,idx_k)*E_N*c_tilde_seg(:,k));
%                     cost = cost -h_N_tilde'*(T_N_tilde(:,:,idx_k)*c_tilde_seg(:,k));
                    cost = cost - 50000*h_N_tilde'*(T_N_tilde(:,:,idx_k)*c_tilde_seg(:,k));
                end
                minimize cost
            
                subject to
                    if i == 1
                        c_seg(:,1) == c0
                        c_tilde_seg(:,1) == c_tilde0
                    else
                        c_seg(:,1) == c_seg_prev(:,end)
                        c_tilde_seg(:,1) == c_tilde_seg_prev(:,end)
                    end
                    
                    for k = 1:num_elem-1
                        idx_k = (num_elem-1)*(i-1)+k;
                        B_N_all\T_N_all(:,:,idx_k)*( E*(c_seg(:,k+1)-c_seg(:,k)) + (tk(idx_k+1)-tk(idx_k))*(M(:,:,idx_k)*c_seg(:,k) - P*kron(eye(prod(N_aug+1)),lambda(:,idx_k))*c_seg(:,k)) ) <= 0
                        B_2N\T_2N(:,:,idx_k+1)*(E_2N*c_tilde_seg(:,k+1) - E_2N_poly*c_seg(:,k+1) - P_2N*kron(eye(prod(N+1)), lambda_tilde(:,idx_k))*c_seg(:,k+1)) <= 0
                    end
            cvx_end
            
            c_seg_prev = c_seg;
            c_tilde_seg_prev = c_tilde_seg;
            
            c(:,(num_elem-1)*(i-1) + (1:num_elem)) = c_seg;
            c_tilde(:,(num_elem-1)*(i-1) + (1:num_elem)) = c_tilde_seg;
            
%             cvx_begin
%                 cvx_precision best
% %                 cvx_quiet true
%             
%                 variable lambda_seg( prod(N_dyn+1), num_elem-1 )
%                 variable lambda_tilde_seg( prod(N+1), num_elem-1 )
%                 
%                 subject to
%                     for k = 1:num_elem-1
%                         idx_k = (num_elem-1)*(i-1)+k;
%                         B_N_all\T_N_all(:,:,idx_k)*( E*(c(:,idx_k+1)-c(:,idx_k)) + (tk(idx_k+1)-tk(idx_k))*M(:,:,idx_k)*c(:,idx_k) - P*kron(c(:,idx_k), eye(prod(N_dyn+1)))*lambda_seg(:,k) ) <= 0
%                         -B_N_dyn\T_N_dyn(:,:,idx_k)*lambda_seg(:,k) <= 0
%                         
%                         B_2N\T_2N(:,:,idx_k+1)*(E_2N*c_tilde(:,idx_k+1) - E_2N_poly*c(:,idx_k+1) - P_2N*kron(c(:,idx_k+1), eye(prod(N+1)))*lambda_tilde_seg(:,k)) <= 0
%                         -B_N\T_N(:,:,idx_k)*lambda_tilde_seg(:,k) <= 0
%                     end
%             cvx_end
%             
%             lambda(:,(num_elem-1)*(i-1) + (1:(num_elem-1))) = lambda_seg;
%             lambda_tilde(:,(num_elem-1)*(i-1) + (1:(num_elem-1))) = lambda_tilde_seg;
        end
        
%         for i = 1:num_seg
%             cvx_begin
%                 cvx_precision best
% %                 cvx_quiet true
%             
%                 variable lambda_seg( prod(N_dyn+1), num_elem-1 )
%                 variable lambda_tilde_seg( prod(N+1), num_elem-1 )
%                 
%                 cost = 0.0;
%                 for k = 1:num_elem-1
%                     idx_k = (num_elem-1)*(i-1)+k;
%                     cost = cost + (E*(c(:,idx_k+1)-c(:,idx_k)) + (tk(idx_k+1)-tk(idx_k))*M(:,:,idx_k)*c(:,idx_k) - P*kron(c(:,idx_k), eye(prod(N_dyn+1)))*lambda_seg(:,k))'*(E*(c(:,idx_k+1)-c(:,idx_k)) + (tk(idx_k+1)-tk(idx_k))*M(:,:,idx_k)*c(:,idx_k) - P*kron(c(:,idx_k), eye(prod(N_dyn+1)))*lambda_seg(:,k)) +...
%                         + (E_2N*c_tilde(:,idx_k+1) - E_2N_poly*c(:,idx_k+1) - P_2N*kron(c(:,idx_k+1), eye(prod(N+1)))*lambda_tilde_seg(:,k))'*(E_2N*c_tilde(:,idx_k+1) - E_2N_poly*c(:,idx_k+1) - P_2N*kron(c(:,idx_k+1), eye(prod(N+1)))*lambda_tilde_seg(:,k));
%                 end
%                 minimize cost
%                 
%                 subject to
%                     for k = 1:num_elem-1
%                         idx_k = (num_elem-1)*(i-1)+k;
%                         B_N_all\T_N_all(:,:,idx_k)*( E*(c(:,idx_k+1)-c(:,idx_k)) + (tk(idx_k+1)-tk(idx_k))*M(:,:,idx_k)*c(:,idx_k) - P*kron(c(:,idx_k), eye(prod(N_dyn+1)))*lambda_seg(:,k) ) <= 0
%                         -B_N_dyn\T_N_dyn(:,:,idx_k)*lambda_seg(:,k) <= 0
%                         
%                         B_2N\T_2N(:,:,idx_k+1)*(E_2N*c_tilde(:,idx_k+1) - E_2N_poly*c(:,idx_k+1) - P_2N*kron(c(:,idx_k+1), eye(prod(N+1)))*lambda_tilde_seg(:,k)) <= 0
%                         -B_N\T_N(:,:,idx_k)*lambda_tilde_seg(:,k) <= 0
%                     end
%             cvx_end
% %             lambda_tilde_seg
%             lambda(:,(num_elem-1)*(i-1) + (1:(num_elem-1))) = lambda_seg;
%             lambda_tilde(:,(num_elem-1)*(i-1) + (1:(num_elem-1))) = lambda_tilde_seg;
%         end

%     cvx_begin
%         variable c_tilde( (sys.Nx)*(sys.Nx+1)/2+1, length(tk) )
%     
%         cost2 = 0.0;
%         for k = 2:length(tk)
%             cost2 = cost2 - h_N_tilde'*(T_N_tilde(:,:,k)*c_tilde(:,k));
%         end
%         minimize cost2
%     
%         subject to
%             c_tilde(:,1) == c_tilde0
%             for k = 1:length(tk)-1
%                 B_N\T_N(:,:,k+1)*(E_N*c_tilde(:,k+1) - c(:,k+1)) <= 0
%             end
%     cvx_end

%         for k = 1:length(tk)-1
%             
%             cvx_begin
%                 cvx_precision best
%                 cvx_quiet true
%             
%                 variable cnext( prod(N+1), 1 )
%             
%                 cost = -100*h_N'*(T_N(:,:,k+1)*cnext);
%                 minimize cost
%             
%                 subject to
%                     cprev = c(:,k);
%                     B\T_N_all(:,:,k)*( E*(cnext-cprev) + (tk(k+1)-tk(k))*(M(:,:,k)*cprev - P*kron(eye(prod(N_aug+1)),lambda_1(:,k))*cprev) ) <= 0
%             cvx_end
%             c(:,k+1) = cnext;
%         end
%     else
%     %% Solve for c
%     cvx_begin
%         variable c( prod(N+1), length(tk) )
%     
%         cost = 0.0;
%         for k = 2:length(tk)
%             cost = cost - 100*h_N'*(T_N(:,:,k)*c(:,k));
%         end
%         minimize cost
%     
%         subject to
%             c(:,1) == c0
%             for k = 1:length(tk)-1
%                 B\T_N_all(:,:,k+1)*( E*(c(:,k+1)-c(:,k)) + (tk(k+1)-tk(k))*(M(:,:,k)*c(:,k) - P*kron(eye(prod(N_aug+1)),lambda_1(:,k))*c(:,k)) ) <= 0
%             end
%     cvx_end
%     end

    %% Solve for lambda
    cvx_begin
        cvx_precision best
        cvx_quiet true
        
        variable lambda( prod(N_dyn+1), length(tk)-1 )
        variable lambda_tilde( prod(N+1), length(tk)-1 )
        
%         cost = 0.0;
%         for k = 1:length(tk)-1
%             cost = cost + (E*(c(:,k+1)-c(:,k))/(tk(k+1)-tk(k)) + M(:,:,k)*c(:,k) - P*kron(c(:,k), eye(prod(N_dyn+1)))*lambda(:,k))'*(E*(c(:,k+1)-c(:,k))/(tk(k+1)-tk(k)) + M(:,:,k)*c(:,k) - P*kron(c(:,k), eye(prod(N_dyn+1)))*lambda(:,k)) +...
%                 (E_2N*c_tilde(:,k+1) - E_2N_poly*c(:,k+1) - P_2N*kron(c(:,k+1), eye(prod(N+1)))*lambda_tilde(:,k))'*(E_2N*c_tilde(:,k+1) - E_2N_poly*c(:,k+1) - P_2N*kron(c(:,k+1), eye(prod(N+1)))*lambda_tilde(:,k));
%         end
%         minimize cost
        
        subject to
            for k = 1:length(tk)-1
                B_N_all\T_N_all(:,:,k)*( E*(c(:,k+1)-c(:,k))/(tk(k+1)-tk(k)) + M(:,:,k)*c(:,k) - P*kron(c(:,k), eye(prod(N_dyn+1)))*lambda(:,k) ) <= 0
%                 -B_N_dyn\T_N_dyn(:,:,k)*lambda(:,k) <= 0
                
                B_2N\T_2N(:,:,k+1)*(E_2N*c_tilde(:,k+1) - E_2N_poly*c(:,k+1) - P_2N*kron(c(:,k+1), eye(prod(N+1)))*lambda_tilde(:,k)) <= 0
%                 -B_N\T_N(:,:,k)*lambda_tilde(:,k) <= 0
            end
    cvx_end
%     lambda
%     lambda_tilde
    %% Post processes
%     c_tilde_conv = c_tilde.*repmat( ones((sys.Nx+1)*(sys.Nx+2)/2,1) - 0.5*quad_idx_skew, [1,length(tk)] );
    c_tilde_conv = c_tilde.*repmat( ones((sys.Nx)*(sys.Nx+1)/2+1,1) - 0.5*quad_idx_skew, [1,length(tk)] );
    r = c_tilde_conv(1,:);
%     d = c_tilde_conv(quad_idx_1st_order, :);
    vechA = c_tilde_conv(quad_idx_2nd_order, :);
    
    % domain to be updated
    region_new = region_;
    for k = 1:length(tk)
        r_ = r(k);
%         d_ = d(:,k);
        A_ = reshape(Utils.duplication_matrix(sys.Nx)*vechA(:,k), [sys.Nx,sys.Nx]);
        
%         q_ = -A_\d_;
%         Q_ = (d_'*(A_\d_)-r_)*inv(A_);
        q_ = zeros(sys.Nx,1);
        Q_ = -r_*inv(A_);
        
        region_new.q(:,k) = q_;
        region_new.Q(:,:,k) = Q_;
    end
    
    % visualization
    if args.plot_cost
        figure(213)
        set(gcf, 'position', [681, 557, 560, 420])
        cla; hold on; grid on;
        for k = round(linspace(1,length(tk),7))
%             tmp = Q_ltv(:,:,k)^(1/2) * Utils.Sphere(1,200).x;
%             plot3(tmp(1,:), tk(k)*ones(1,size(tmp,2)), tmp(2,:), 'g', 'linewidth', 2)
            
            if sum(strcmp(fieldnames(args), 'frs')) == 1
                tmp = args.frs{k};
                plot3(tmp(1,:), tk(k)*ones(1,size(tmp,2)), tmp(2,:), 'k*', 'linewidth', 2);
            end
            
            Vp = Polynomial(c(:,k), Polynomial.integer_grid(N));
            tmp = Utils.get_level_set( args.grid, Vp.eval(args.grid), 0.0 );
            plot3(tmp(1,:), tk(k)*ones(1,size(tmp,2)), tmp(2,:), 'r', 'linewidth', 2)
            
            tmp = region_new.Q(:,:,k)^(1/2) * Utils.Sphere(1,200).x;
            plot3(tmp(1,:), tk(k)*ones(1,size(tmp,2)), tmp(2,:), 'b--', 'linewidth', 2)
        end
        xlabel('$x_1$')
        zlabel('$x_2$')
        ylabel('$t$ [s]')
        view([122,19])
        drawnow
        
        
        vol = 0.0;
        for k = 2:length(tk)
%             Vp = Polynomial(c(:,k), Polynomial.integer_grid(N));
%             tmp = Utils.get_level_set( args.grid, Vp.eval(args.grid), 0.0 );
%             vol = vol + log( area(polyshape(tmp')) );
            vol = vol + log(det( region_new.Q(:,:,k) ));
        end
        dec_rate = (cost_hist(end) - vol) / cost_hist(end);
        disp(['Iteration #', num2str(iter),...
            ': cost = ', num2str(vol),...
            ', dcost = ', num2str(cost_hist(end)-vol),...
            ', rate = ', num2str(dec_rate)]);
        
        cost_hist = [cost_hist, vol];
        rate_hist = [rate_hist, dec_rate];
        
        
        figure(123)
        set(gcf, 'position', [1242, 557, 560, 420])
        subplot(2,1,1)
        title('$\textbf{Convergence characteristics of SLP}$')
        cla; hold on; grid on;
        plot(cost_hist,'b*-')
        axis tight;
        ylabel('$cost$')
        
        subplot(2,1,2)
        cla; hold on; grid on;
        plot(rate_hist,'b*-')
        axis tight;
        ylabel('$rate$')
        xlabel('Iterations')
        drawnow
    end
    
    %% Convergence check & domain update
    if abs(dec_rate) < args.ftol
        disp(['Converged (rate = ',num2str(dec_rate),' / ftol = ',num2str(args.ftol),')'])
        break;
    end
    region(iter+1) = region_new;

end
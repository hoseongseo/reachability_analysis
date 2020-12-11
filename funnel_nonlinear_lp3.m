function [region, cost_hist, rate_hist, c, s] = funnel_nonlinear_lp3(sys, t, Q0, N, N_lambda, N_gamma, args)
% sys: dynamics
% tk: time horizon of interest (discretized)
% Q0: shape matrix for the initial set
% N: target order of the approximated value function
% N_lambda: multiplier polynomial for polynomial approximation
% N_gamma: multiplier polynomial for quadratic approximation
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
Q_ltv = funnel_ltv(sys, t, Q0);
region = struct('q', zeros(sys.Nx,length(t)), 'Q', Q_ltv);

%%%%%% degrees of the approximated value functions
N_a = [N, zeros(1,sys.Nw)]; % augmentation
fp = sys.f(xp, zeros(sys.Nu,1), wp, t(1));
N_f = max(reshape([fp.order], [numel([fp.order])/(sys.Nx+sys.Nw), (sys.Nx+sys.Nw)]), [], 1);

N_max = max(N_a + N_f, N_a + N_lambda); % total order of polynmial approximation 
N_max_quad = N + N_gamma; % total order of quadratic approximation


poly_degree = Polynomial.integer_grid(N);
quad_degree = poly_degree((sum(poly_degree, 2) == 2) | (sum(poly_degree, 2) == 0), :);
quad_idx_2nd_order = sum(quad_degree,2) == 2;
quad_idx_skew = max(quad_degree,[],2) == 1;


E_N = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a));
E_N_N_f = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a + N_f));
E_N_N_lambda = Polynomial.expand_matrix(N_max, Polynomial.integer_grid(N_a + N_lambda));
E_N_max_quad = Polynomial.expand_matrix(N_max_quad, quad_degree);
E_N_max_quad_poly = Polynomial.expand_matrix(N_max_quad, poly_degree);

B_N_max = Polynomial.bernstein_transform_matrix(N_max);
B_N_max_quad = Polynomial.bernstein_transform_matrix(N_max_quad);
B_lambda = Polynomial.bernstein_transform_matrix(N_lambda);
B_gamma = Polynomial.bernstein_transform_matrix(N_gamma);

P_N_N_f = Polynomial.product_matrix(N_a, N_f);
P_N_N_lambda = Polynomial.product_matrix(N_a, N_lambda);
P_N_N_gamma = Polynomial.product_matrix(N, N_gamma);


D = Polynomial.derivative_matrix(N);

% initial coefficients
V0 = xp'*inv(Q0)*xp - 1;
c0 = Polynomial.expand_matrix(N, V0.order(:,1:sys.Nx))*V0.coeff;
c_quad0 = c0((sum(poly_degree, 2) == 2) | (sum(poly_degree, 2) == 0));

M = zeros( prod(N_a+N_f+1), prod(N+1), length(t)-1 );
for k = 1:length(t)-1
    mat = zeros(prod(N_a+N_f+1), prod(N+1));
    f_ = sys.f(xp, zeros(sys.Nu,1), wp, t(k));
    for j = 1:sys.Nx
        mat = mat + P_N_N_f * kron( eye(prod(N+1)), Polynomial.expand_matrix(N_f, f_(j).order)*f_(j).coeff ) * D(:,:,j);
    end
    M(:,:,k) = mat;
end


lambda = zeros( prod(N_lambda+1), length(t)-1 ); % space for lambda
gamma = zeros( prod(N_gamma+1), length(t)-1 ); 

% volume of the initial guess
vol = 0.0;
for k = 2:length(t)
    vol = vol + log(det( region.Q(:,:,k) ));
end
cost_hist = vol;
rate_hist = [];


T_N_max = zeros( prod(N_max+1), prod(N_max+1), length(t) );
T_N_max_quad = zeros( prod(N_max_quad+1), prod(N_max_quad+1), length(t) );
T_lambda = zeros( prod(N_lambda+1), prod(N_lambda+1), length(t));
T_gamma = zeros( prod(N_gamma+1), prod(N_gamma+1), length(t));

[~, full_order_quad] = Polynomial.affine_transform_matrix2(quad_degree, Q0^(1/2), zeros(sys.Nx,1));
T_quad = zeros( size(full_order_quad,1), size(quad_degree,1), length(t) );
h_quad = prod((1./(full_order_quad+1)).*( ones(size(full_order_quad)).^(full_order_quad+1) - (-ones(size(full_order_quad))).^(full_order_quad+1) ),2);


region_ = region(end);
for k = 1:length(t)
    sqrtQ_ = region_.Q(:,:,k)^(1/2);
    a_ = sqrt(sum(sqrtQ_.^2,2));
    b_ = region_.q(:,k);
    
    T_N_max(:,:,k) = Polynomial.domain_transform_matrix( N_max, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
    T_quad(:,:,k) = Polynomial.affine_transform_matrix2( quad_degree, sqrtQ_, b_ );
    T_N_max_quad(:,:,k) = Polynomial.domain_transform_matrix( N_max_quad, a_, b_ );
    T_lambda(:,:,k) = Polynomial.domain_transform_matrix( N_lambda, [a_; ones(sys.Nw,1)], [b_;zeros(sys.Nw,1)] );
    T_gamma(:,:,k) = Polynomial.domain_transform_matrix( N_gamma, a_, b_ );
end


% MAIN LOOP
num_elem = 11; % number of elements per segment
num_seg = (length(t)-1)/(num_elem-1);

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
        c = zeros(prod(N+1), length(t));
        c_quad = zeros((sys.Nx)*(sys.Nx+1)/2+1, length(t));
                
        for i = 1:num_seg
            cvx_begin
                cvx_precision default
%                 cvx_quiet true
            
                variable c_seg( prod(N+1), num_elem )
                variable c_quad_seg( (sys.Nx)*(sys.Nx+1)/2+1, num_elem )
                
                cost = 0.0;
                for k = 2:num_elem
                    idx_k = (num_elem-1)*(i-1)+k;
                    cost = cost - h_quad'*(T_quad(:,:,idx_k)*c_quad_seg(:,k));
                end
                minimize cost
            
                subject to
                    if i == 1
                        c_seg(:,1) == c0
                        c_quad_seg(:,1) == c_quad0
                    else
                        c_seg(:,1) == c_seg_prev(:,end)
                        c_quad_seg(:,1) == c_quad_seg_prev(:,end)
                    end
                    
                    for k = 1:num_elem-1
                        idx_k = (num_elem-1)*(i-1)+k;
                        B_N_max\T_N_max(:,:,idx_k)*( E_N*(c_seg(:,k+1)-c_seg(:,k)) + (t(idx_k+1)-t(idx_k))*(E_N_N_f*M(:,:,idx_k)*c_seg(:,k) - E_N_N_lambda*P_N_N_lambda*kron(eye(prod(N+1)),lambda(:,idx_k))*c_seg(:,k)) ) <= 0
                        
                        B_N_max_quad\T_N_max_quad(:,:,idx_k+1)*(E_N_max_quad*c_quad_seg(:,k+1) - E_N_max_quad_poly*c_seg(:,k+1) - P_N_N_gamma*kron(eye(prod(N+1)), gamma(:,idx_k))*c_seg(:,k+1)) <= 0
                    end
            cvx_end
            
            c_seg_prev = c_seg;
            c_quad_seg_prev = c_quad_seg;
            
            c(:,(num_elem-1)*(i-1) + (1:num_elem)) = c_seg;
            c_quad(:,(num_elem-1)*(i-1) + (1:num_elem)) = c_quad_seg;
            
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
        cvx_precision default
%         cvx_quiet true
        
        variable lambda( prod(N_lambda+1), length(t)-1 )
        variable gamma( prod(N_gamma+1), length(t)-1 )
        variable c_quad_tmp( (sys.Nx)*(sys.Nx+1)/2+1, length(t) ) 
        
        cost = 0.0;
        for k = 1:length(t)-1
            cost = cost - h_quad'*(T_quad(:,:,k+1)*c_quad_tmp(:,k+1));
        end
        minimize cost
        
        subject to
%             c_quad_tmp(:,1) == c_quad0
            
            for k = 1:length(t)-1
                B_N_max\T_N_max(:,:,k)*( E_N*(c(:,k+1)-c(:,k)) + (t(k+1)-t(k))*(E_N_N_f*M(:,:,k)*c(:,k) - E_N_N_lambda*P_N_N_lambda*kron(c(:,k), eye(prod(N_lambda+1)))*lambda(:,k)) ) <= 0
                -B_lambda\T_lambda(:,:,k)*lambda(:,k) <= 0
                
%                 B_N_max_quad\T_N_max_quad(:,:,k+1)*(E_N_max_quad*c_quad(:,k+1) - E_N_max_quad_poly*c(:,k+1) - P_N_N_gamma*kron(c(:,k),eye(prod(N_gamma+1)))*gamma(:,k)) <= 0
                B_N_max_quad\T_N_max_quad(:,:,k+1)*(E_N_max_quad*c_quad_tmp(:,k+1) - E_N_max_quad_poly*c(:,k+1) - P_N_N_gamma*kron(c(:,k),eye(prod(N_gamma+1)))*gamma(:,k)) <= 0
                -B_gamma\T_gamma(:,:,k)*gamma(:,k) <= 0
            end
    cvx_end
%     lambda
%     lambda_tilde
    %% Post processes
    c_quad_conv = c_quad.*repmat( ones((sys.Nx)*(sys.Nx+1)/2+1,1) - 0.5*quad_idx_skew, [1,length(t)] );
    r = c_quad_conv(1,:);
    vechA = c_quad_conv(quad_idx_2nd_order, :);
    
    % domain to be updated
    region_new = region_;
    for k = 1:length(t)
        r_ = r(k);
        A_ = reshape(Utils.duplication_matrix(sys.Nx)*vechA(:,k), [sys.Nx,sys.Nx]);
        
        q_ = zeros(sys.Nx,1);
        Q_ = -r_*inv(A_);
        
        region_new.q(:,k) = q_;
        region_new.Q(:,:,k) = Q_;
    end
    
    c_quad_conv_tmp = c_quad_tmp.*repmat( ones((sys.Nx)*(sys.Nx+1)/2+1,1) - 0.5*quad_idx_skew, [1,length(t)] );
    r_tmp = c_quad_conv_tmp(1,:);
    vechA_tmp = c_quad_conv_tmp(quad_idx_2nd_order, :);
    
    % domain to be updated
    region_new_tmp = region_;
    for k = 1:length(t)
        r_ = r_tmp(k);
        A_ = reshape(Utils.duplication_matrix(sys.Nx)*vechA_tmp(:,k), [sys.Nx,sys.Nx]);
        
        q_ = zeros(sys.Nx,1);
        Q_ = -r_*inv(A_);
        
        region_new_tmp.q(:,k) = q_;
        region_new_tmp.Q(:,:,k) = Q_;
    end
    
    % visualization
    if args.plot_cost
        figure(213)
        set(gcf, 'position', [681, 557, 560, 420])
        cla; hold on; grid on;
        for k = round(linspace(1,length(t),5))
            tmp = Q_ltv(:,:,k)^(1/2) * Utils.Sphere(1,200).x;
            plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'g', 'linewidth', 2)
            
            if sum(strcmp(fieldnames(args), 'frs')) == 1
                tmp = args.frs{k};
                plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'k*', 'linewidth', 2);
            end
            
            Vp = Polynomial(c(:,k), Polynomial.integer_grid(N));
            tmp = Utils.get_level_set( args.grid, Vp.eval(args.grid), 0.0 );
            plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'r', 'linewidth', 2)
            
            tmp = region_new.Q(:,:,k)^(1/2) * Utils.Sphere(1,200).x;
            plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'b--', 'linewidth', 2)
            
            tmp = region_new_tmp.Q(:,:,k)^(1/2) * Utils.Sphere(1,200).x;
            plot3(tmp(1,:), t(k)*ones(1,size(tmp,2)), tmp(2,:), 'm.', 'linewidth', 2)
        end
        xlabel('$x_1$')
        zlabel('$x_2$')
        ylabel('$t$ [s]')
        view([122,19])
        drawnow
        
        
        vol = 0.0;
        for k = 2:length(t)
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
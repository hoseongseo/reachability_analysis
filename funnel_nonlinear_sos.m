function [region, cost_hist, rate_hist, ctime_hist] = funnel_nonlinear_sos( sys, tk, Q0, args )
% sys: dynamics
% tk: time horizon of interest (discretized)
% Q0: shape matrix for the initial set
% args: extra arguments detailed as
if nargin < 4
    args.max_iter = 1;
    args.rho = 3.0;
    args.ftol = 5e-4;
    args.plot_cost = false;
end

solver_opt.solver = 'sedumi';

nx = sys.Nx;
nw = sys.Nw;

N = length(tk);

x = mpvar('x', [nx,1]);
w = mpvar('w', [nw,1]);
vars = [x;w];

% state bound
gx = 1 - x'*inv(Q0)*x;

% disturbance bound
gw = cell(1,nw);
for j = 1:nw
    gw{j} = 1 - (w(j)*w(j));
end

% initial guess on level
rho_init = exp(args.rho*tk);

% preparation
step2_solved = false;
cost_hist = [];
rate_hist = [];
ctime_hist = [];

order = [0,1,2];
res = 0;

region = struct('step1', [], 'step2', []);
Q_step1 = zeros(nx,nx,N);
Q_step2 = zeros(nx,nx,N);

% MAIN LOOP
for iter = 1:args.max_iter
    %% STEP 1 : Lagrange multiplier polynomial
    step1 = sosprogram(vars);
    
    % Variable 1 : Lyapunov polynomial for initial condition
    [step1,L0] = sospolyvar(step1,monomials(x,order));
    
    % Varables 2 : Lyapunov polynomials for the value function
    L = cell(1,N);
    for k = 1:N
        [step1,L{k}] = sospolyvar(step1,monomials(vars,order));
    end
    
    % Varables 3 : Lyapunov polynomials for outer ellipsoids
    Le = cell(1,N);
    for k = 1:N
        [step1,Le{k}] = sospolyvar(step1,monomials(x,order));
    end
    
    % Variables 4 : Shape matrix for outer ellipsoid
    S = cell(1,N);
    for k = 1:N
        [step1, S{k}] = sospolymatrixvar(step1,monomials(vars,0),[nx,nx],'symmetric');
    end
    
    % Varables 5 : Lyapunov polynomials for disturbance
    Lw_1 = cell(nw,N);
    for k = 1:N
        for j = 1:nw
            [step1,Lw_1{j,k}] = sospolyvar(step1,monomials(vars,order));
        end
    end
    
    % Constraint 1 : Initial condition
    % rho(t_1) - V_1(x) - L0(x)*g0(x) is SOS, where g0(x) = 1 - x'*inv(Q0)*x
    if step2_solved
        V0 = sosgetsol(step2, V{1});
        rho0 = sosgetsol(step2, rho{1});
    else
        V0 = x'*inv(Q0)*x;
        rho0 = rho_init(1);
    end
    step1 = sosineq(step1, rho0 - V0 - L0*gx);
    step1 = sosineq(step1, L0); % L0 should be SOS itself
    
    obj1 = polynomial(0); % objective function
    for k = 1:N
        if step2_solved
            rho_ = sosgetsol(step2, rho{k});
            V_ = sosgetsol(step2, V{k});
        else
            rho_ = rho_init(k);
            V_ = x'*inv(Q0)*x;
        end
        
        S_ = S{k};
        
        % Constraints 2 : Ellipsoid containtment constraint
        % 1 - x'*S_k*x - Le_k*(p(t_k) - V_k(x)) is SOS for all k = 1:N
        step1 = sosineq(step1, 1 - x'*S_*x - Le{k}*(rho_ - V_));
        
        % Constraints 3 : Ensuring ellipsoid is SOS
        % S_k is SOS
        step1 = sosmatrixineq(step1, S_);
        
        % Constraints 4 : Le is SOS
        step1 = sosineq(step1, Le{k});
        
        % Constraints 5 : Value gradients
        if k > 1
            dt_ = tk(k) - tk(k-1);
            drho_ = (rho_ - rhoprev_) / dt_;
            dV_dt_ = (V_ - Vprev_) / dt_; % partial derivative of V wrt time
            
            if step2_solved
                dV_dx_ = [];
                for i = 1:nx
                    dV_dx_ = [dV_dx_; diff(V_,x(i))];
                end
            else
                dV_dx_ = 2*inv(Q0)*x;
            end
            dV_ = dV_dt_ + dV_dx_'*sys.f(x, zeros(sys.Nu,1), w, tk(k));
            
            
%             expr = drho_ - dV_ - L{k}*(V_ - rho_);
            expr = drho_ - dV_ - L{k}*(rho_ - V_);
            for j = 1:nw
                expr = expr - Lw_1{j,k}*gw{j};
            end
            step1 = sosineq(step1, expr);
        end
        
        % Constraint 6 : Lw is SOS
        for j = 1:nw
            step1 = sosineq(step1, Lw_1{j,k});
        end
        
        % update cost function
        if step2_solved
            obj1 = obj1 - trace(S_*inv(S_prev(:,:,k)));
        else
            obj1 = obj1 - trace(S_*(Q0*rho_init(k)));
        end
%         obj1 = obj1 - trace(S_);
        Vprev_ = V_;
        rhoprev_ = rho_;
    end
    
    step1 = sossetobj(step1,obj1);
    
%     tic;
    step1 = sossolve(step1,solver_opt);
%     ctime1 = toc;

    %% STEP 2 : Value and level (L and Lem fixed)
    step2 = sosprogram(vars);
    
    % Variable 1 : Lyapunov polynomial for initial condition
    [step2,L1] = sospolyvar(step2,monomials(x,order)); % should this be also quadratic??
    
    % Varables 2 : Value polynomial
    V = cell(1,N);
    for k = 1:N
        [step2,V{k}] = sospolyvar(step2,monomials(x,order));
    end
    
    % Varables 3 : Level polynomial
    rho = cell(1,N);
    for k = 1:N
        [step2,rho{k}] = sospolyvar(step2,monomials(vars,0));
    end
    
    % Variables 4 : Shape matrix for outer ellipsoid (same with S in STEP 1)
    R = cell(1,N);
    for k = 1:N
        [step2, R{k}] = sospolymatrixvar(step2,monomials(vars,0),[nx,nx],'symmetric');
    end
    
    % Varables 5 : Lyapunov polynomials for disturbance
    Lw_2 = cell(nw,N);
    for k = 1:N
        for j = 1:nw
            [step2,Lw_2{j,k}] = sospolyvar(step2,monomials(vars,order));
        end
    end
    
    % Constraint 1 : Initial condition
    % rho(t_1) - V_1(x) - L0(x)*g0(x) is SOS, where g0(x) = 1 - x'*inv(Q0)*x
    step2 = sosineq(step2, rho{1} - V{1} - L1*gx);
    step2 = sosineq(step2, L1); % L0 should be SOS itself !
    
    obj2 = polynomial(0); % objective function
    for k = 1:N
        rho_ = rho{k};
        V_ = V{k};
        S_ = R{k};
        
        % Constraints 2 : Ellipsoid containtment constraint
        % 1 - x'*S_k*x - Le_k*(p(t_k) - V_k(x)) is SOS for all k = 1:N
        step2 = sosineq(step2, 1 - x'*S_*x - sosgetsol(step1, Le{k})*(rho_ - V_));
        
        % Constraints 3 : Ensuring ellipsoid is SOS
        step2 = sosmatrixineq(step2, S_);
        
        % Constraints 4 : Value gradients
        if k > 1
            dt_ = tk(k) - tk(k-1);
            drho_ = (rho_ - rhoprev_) / dt_;
            dV_dt_ = (V_ - Vprev_) / dt_; % partial derivative of V wrt time
            dV_dx_ = [];
            for i = 1:nx
                dV_dx_ = [dV_dx_; diff(V_,x(i))];
            end
            dV_ = dV_dt_ + dV_dx_'*sys.f(x, zeros(sys.Nu,1), w, tk(k));
            
%             expr = drho_ - dV_ - sosgetsol(step1, L{k})*(V_ - rho_);
            expr = drho_ - dV_ - sosgetsol(step1, L{k})*(rho_ - V_);
            for j = 1:nw
                expr = expr - Lw_2{j,k}*gw{j};
            end
            step2 = sosineq(step2, expr);
        end
        
        % Constraint 5 : Lw is SOS
        for j = 1:nw
            step2 = sosineq(step2, Lw_2{j,k});
        end
        
        % update cost function
        if step2_solved
            obj2 = obj2 - trace(S_*inv(R_prev(:,:,k)));
        else
            obj2 = obj2 - trace(S_*(Q0*rho_init(k)));
        end
%         obj2 = obj2 - trace(S_);
        Vprev_ = V_;
        rhoprev_ = rho_;
    end
    
    step2 = sossetobj(step2,obj2);
%     tic
    step2 = sossolve(step2,solver_opt);
%     ctime2 = toc;
    
    %% check cost and dcost
    cost = 0;
    for k = 1:N
        Q_step1(:,:,k) = inv(double(sosgetsol(step1, S{k})));
        Q_step2(:,:,k) = inv(double(sosgetsol(step2, R{k})));
        cost = cost + log(det( Q_step2(:,:,k) ));
    end
    region_ = struct('step1', Q_step1, 'step2', Q_step2);
    
    if ~step2_solved
        disp(['Iteration #', num2str(iter),...
            ': cost = ', num2str(cost)])
        cost_hist = [cost_hist, cost];
%         ctime_hist = [ctime_hist, [ctime1; ctime2]];
    else
        dec_rate = (cost_prev - cost) / cost_prev;
        disp(['Iteration #', num2str(iter),...
            ': cost = ', num2str(cost),...
            ', dcost = ', num2str(cost_prev-cost),...
            ', rate = ', num2str(dec_rate)])
        rate_hist = [rate_hist, dec_rate];
        cost_hist = [cost_hist, cost];
%         ctime_hist = [ctime_hist, [ctime1; ctime2]];
            
        if abs(dec_rate) < args.ftol 
            % converged, accept the last result
            res = 1; 
            disp(['Converged (rate = ',num2str(dec_rate),' / ftol = ',num2str(args.ftol),')'])
            break;
        end
%         if dec_rate > 0
%             % cost increased, reject the last result
%             res = 2;
%             disp('Terminate due to increase of cost')
%             break;
%         end
    end
    region(iter) = region_;
    
    % save the result of the current step as S_prev
    S_prev = zeros(nx,nx,N);
    for k = 1:N
        S_prev(:,:,k) = double(sosgetsol(step1, S{k}));
    end
    
    % save the result of the current step as R_prev
    R_prev = zeros(nx,nx,N);
    for k = 1:N
        R_prev(:,:,k) = double(sosgetsol(step2, R{k}));
    end
    
    cost_prev = cost;
    step2_solved = true;
    
    if args.plot_cost
        figure(124)
        subplot(2,1,1)
        title('$\textbf{Convergence characteristics of SOS program}$')
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
        
%         Sp = Math.Sphere(1,200);
%         
%         figure(2213)
%         cla; hold on; grid on;
%         for k = 1:length(sys.tN)
%             tmp = sys.xN(1:2,k) + Q_step1(1:2,1:2,k)^(1/2) * Sp.x;
%             plot(tmp(1,:), tmp(2,:), 'b', 'linewidth', 2)
%             
%             tmp = sys.xN(1:2,k) + Q_step2(1:2,1:2,k)^(1/2) * Sp.x;
%             plot(tmp(1,:), tmp(2,:), 'r', 'linewidth', 2)
%         end
% %         for k = 1:size(args.xN,2)
% %             tmp = args.xN(1:2,k) + args.xChar(1:2,:,k);
% %             plot(tmp(1,:), tmp(2,:), 'k.', 'linewidth', 2)
% %         end
%         drawnow
    end
end

if iter == args.max_iter
    res = 3; % maximum iteration
    disp('Maximum iteration reached')
end

if (res == 1) || (res == 3)
    region(iter) = region_;
end
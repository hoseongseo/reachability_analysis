function Q_funnel = invariant_funnel( t, xNominal, f_poly, Q, wMax, rho, maxIter, ftol )
solver_opt.solver = 'sedumi';

nx = size(Q,1);
nw = length(wMax);

N = length(t);

x = mpvar('x', [nx,1]);
w = mpvar('w', [nw,1]);
vars = [x;w];

% state bound
gx = 1 - x'*inv(Q)*x;

% disturbance bound
gw = cell(1,nw);
for j = 1:nw
    gw{j} = 1 - (w(j)*w(j))/(wMax(j)^2);
end

% initial guess on level
rhoInit = exp(rho*t);

% preparation
step2_solved = false;
costHist = [];
rateHist = [];

order = 2;
for iter = 1:maxIter
    %% STEP 1 : Lagrange multiplier polynomial
    costHist1 = [];
    rateHist1 = [];
    for innerIter = 1:1
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
            V0 = x'*inv(Q)*x;
            rho0 = rhoInit(1);
        end
        step1 = sosineq(step1, rho0 - V0 - L0*gx);
        step1 = sosineq(step1, L0); % L0 should be SOS itself
        
        obj1 = polynomial(0); % objective function
        for k = 1:N
            xNominal_ = xNominal(:,k);
            
            if step2_solved
                rho_ = sosgetsol(step2, rho{k});
                V_ = sosgetsol(step2, V{k});
            else
                rho_ = rhoInit(k);
                V_ = x'*inv(Q)*x;
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
                dt_ = t(k) - t(k-1);
                drho_ = (rho_ - rhoprev_) / dt_;
                dV_dt_ = (V_ - Vprev_) / dt_; % partial derivative of V wrt time
                
                if step2_solved
                    dV_dx_ = [];
                    for i = 1:nx
                        dV_dx_ = [dV_dx_; diff(V_,x(i))];
                    end
                else
                    dV_dx_ = 2*inv(Q)*x;
                end
                dV_ = dV_dt_ + dV_dx_'*f_poly(x,w,xNominal_);
                
                
                expr = drho_ - dV_ - L{k}*(V_ - rho_);
                for j = 1:nw
                    expr = expr - Lw_1{j,k}*gw{j};
                end
                step1 = sosineq(step1, expr);
            end
            
            % Constraint 6 : Lw is SOS
            for j = 1:nw
                step1 = sosineq(step1, Lw_1{j,k});
            end
            
%             % update cost function
%             obj1 = obj1 - trace(S_); %%%%%%% REVISE THIS!!
            if step2_solved
                % obj1 = obj1 - trace(S_*inv(double(sosgetsol(step2, R{k}))));
                obj1 = obj1 - trace(S_*inv(S_prev(:,:,k)));
            else
                obj1 = obj1 - trace(S_*(Q*rhoInit(k)));
            end
%             if innerIter == 1
%                 if step2_solved
%                     obj1 = obj1 - trace(S_*inv(double(sosgetsol(step2, R{k}))));
%                 else
%                     obj1 = obj1 - trace(S_*(Q*rhoInit(k)));
%                 end
%             else
%                 obj1 = obj1 - trace(S_*inv(S_prev(:,:,k)));
%             end
%                 
            Vprev_ = V_;
            rhoprev_ = rho_;
        end
        
        step1 = sossetobj(step1,obj1);
        step1 = sossolve(step1,solver_opt);
        
        cost1 = 0;
        for k = 1:N
            cost1 = cost1 - double(det( sosgetsol(step1, S{k}) ));
        end
        
        figure(444)
        cla; hold on; grid on;
        for k = round(linspace(1,N,11))
            % STEP1
            Q_sos_ = inv(double(sosgetsol(step1, S{k})));
            tmp = xNominal(:,k) + Q_sos_^(1/2)*Math.Sphere(nx-1,100).x;
            plot(tmp(1,:), tmp(2,:), 'b', 'linewidth', 2)
            
            if innerIter > 1
                Q_sos_ = inv(S_prev(:,:,k));
                tmp = xNominal(:,k) + Q_sos_^(1/2)*Math.Sphere(nx-1,100).x;
                plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2)
            end
        end
        
        costHist1 = [costHist1, cost1];
        if innerIter > 1
            dec_rate1 = (cost1_prev - cost1) / cost1_prev;
            rateHist1 = [rateHist1, dec_rate1];
            if abs(dec_rate1) < 1e-3
                disp('Inner loop #1 converged')
                break;
            else
                S_prev = zeros(nx,nx,N);
                for k = 1:N
                    S_prev(:,:,k) = double(sosgetsol(step1, S{k}));
                end
                cost1_prev = cost1;
            end            
        else
            S_prev = zeros(nx,nx,N);
            for k = 1:N
                S_prev(:,:,k) = double(sosgetsol(step1, S{k}));
            end
            cost1_prev = cost1;
        end
%         S_prev = zeros(nx,nx,N);
%         for k = 1:N
%             S_prev(:,:,k) = double(sosgetsol(step1, S{k}));
%         end


        figure(333)
        subplot(2,1,1)
        cla; hold on; grid on;
        plot(costHist1, 'bx-')
        axis tight;
        subplot(2,1,2)
        cla; hold on; grid on;
        plot(rateHist1, 'bx-')
        axis tight;
        drawnow
    end
    %% STEP 2 : Value and level (L and Lem fixed)
    costHist2 = [];
    rateHist2 = [];
    for innerIter = 1:1
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
            xNominal_ = xNominal(:,k);
            
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
                dt_ = t(k) - t(k-1);
                drho_ = (rho_ - rhoprev_) / dt_;
                dV_dt_ = (V_ - Vprev_) / dt_; % partial derivative of V wrt time
                dV_dx_ = [];
                for i = 1:nx
                    dV_dx_ = [dV_dx_; diff(V_,x(i))];
                end
                dV_ = dV_dt_ + dV_dx_'*f_poly(x,w,xNominal_);
                
                expr = drho_ - dV_ - sosgetsol(step1, L{k})*(V_ - rho_);
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
%             obj2 = obj2 - trace(S_);
%             obj2 = obj2 - trace(S_*inv(double(sosgetsol(step1, S{k}))));
            if step2_solved
                obj2 = obj2 - trace(S_*inv(R_prev(:,:,k)));
            else
                obj2 = obj2 - trace(S_*(Q*rhoInit(k)));
            end
%             if innerIter == 1
% %                 obj2 = obj2 - trace(S_*inv(double(sosgetsol(step1, S{k}))));
%                 obj2 = obj2 - trace(S_);
%             else
%                 obj2 = obj2 - trace(S_*inv(R_prev(:,:,k)));
%             end
            Vprev_ = V_;
            rhoprev_ = rho_;
        end
        
        step2 = sossetobj(step2,obj2);
        step2 = sossolve(step2,solver_opt);
        
        cost2 = 0;
        for k = 1:N
            cost2 = cost2 - double(det( sosgetsol(step2, R{k}) ));
        end
        
        figure(4444)
        cla; hold on; grid on;
        for k = round(linspace(1,N,11))
            % STEP1
            Q_sos_ = inv(double(sosgetsol(step2, R{k})));
            tmp = xNominal(:,k) + Q_sos_^(1/2)*Math.Sphere(nx-1,100).x;
            plot(tmp(1,:), tmp(2,:), 'b', 'linewidth', 2)
            
            if innerIter > 1
                Q_sos_ = inv(R_prev(:,:,k));
                tmp = xNominal(:,k) + Q_sos_^(1/2)*Math.Sphere(nx-1,100).x;
                plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2)
            end
        end
        
        costHist2 = [costHist2, cost2];
        if innerIter > 1
            dec_rate2 = (cost2_prev - cost2) / cost2_prev;
            rateHist2 = [rateHist2, dec_rate2];
            if abs(dec_rate2) < 1e-3
                disp('Inner loop #1 converged')
                break;
            else
                R_prev = zeros(nx,nx,N);
                for k = 1:N
                    R_prev(:,:,k) = double(sosgetsol(step2, R{k}));
                end
                cost2_prev = cost2;
            end            
        else
            R_prev = zeros(nx,nx,N);
            for k = 1:N
                R_prev(:,:,k) = double(sosgetsol(step2, R{k}));
            end
            cost2_prev = cost2;
        end

        figure(3333)
        subplot(2,1,1)
        cla; hold on; grid on;
        plot(costHist2, 'bx-')
        axis tight;
        subplot(2,1,2)
        cla; hold on; grid on;
        plot(rateHist2, 'bx-')
        axis tight;
        drawnow
        
%         R_prev = zeros(nx,nx,N);
%         for k = 1:N
%             R_prev(:,:,k) = double(sosgetsol(step2, R{k}));
%         end
    end
    
    cost = 0;
    for k = 1:N
        cost = cost - double(det( sosgetsol(step2, R{k}) ));
    end

    if ~step2_solved
        disp(['Iteration #', num2str(iter),...
            ': cost = ', num2str(cost)])
    else
        dec_rate = (cost_prev - cost) / cost_prev;
        disp(['Iteration #', num2str(iter),...
            ': cost = ', num2str(cost),...
            ', dcost = ', num2str(cost_prev-cost),...
            ', rate = ', num2str(dec_rate)])
        rateHist = [rateHist, dec_rate];
        
        if abs(dec_rate) < ftol
            disp(['Converged (rate = ',num2str(dec_rate),' / ftol = ',num2str(ftol),')'])
            break;
        end
    end
    costHist = [costHist, cost];

    cost_prev = cost;
    step2_solved = true;
    
    %% confirm
    figure(33)
    cla; hold on; grid on;
    for k = round(linspace(1,N,11))
        % STEP1
        Q_sos_ = inv(double(sosgetsol(step1, S{k})));
        tmp = xNominal(:,k) + Q_sos_^(1/2)*Math.Sphere(nx-1,100).x;
        plot(tmp(1,:), tmp(2,:), 'b', 'linewidth', 2)
  
        % STEP2
        Q_sos_ = inv(double(sosgetsol(step2, R{k})));
        tmp = xNominal(:,k) + Q_sos_^(1/2)*Math.Sphere(nx-1,100).x;
        plot(tmp(1,:), tmp(2,:), 'r--', 'linewidth', 2)
    end

    figure(34)
    subplot(2,1,1)
    cla; hold on; grid on;
    plot(costHist,'b*-')
    axis tight;
    ylabel('Cost')
    subplot(2,1,2)
    cla; hold on; grid on;
    plot(rateHist,'b*-')
    axis tight;
    ylabel('Rate')
    xlabel('Iteration')
    drawnow
end

if iter == maxIter
    disp('Maximum iteration reached')
end

Q_funnel = zeros(nx,nx,N);
for k = 1:N
    Q_funnel(:,:,k) = inv(double(sosgetsol(step2, R{k})));
end
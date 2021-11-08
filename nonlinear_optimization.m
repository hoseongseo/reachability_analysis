function Qres = nonlinear_optimization(sys, t, Q0, wMax)

dt = t(2) - t(1);

% discretized dynamics
f = @(x,w)...
    [x(1,:) + dt*(w+3).*x(1,:).*(1-x(2,:));...
    x(2,:) + dt*(w+3).*x(2,:).*(x(1,:)-1)];

qres = zeros(2,length(t));
qres(:,1) = [1.2; 1.1];
Qres = zeros(sys.Nx,sys.Nx,length(t));
Qres(:,:,1) = Q0;

options = optimoptions('fmincon',...
    'Algorithm', 'sqp',...
    'Display', 'none',...
    'MaxFunctionEvaluations', 10000000,...
    'MaxIterations', 10000000);

for k = 1:length(t)-1
    q = qres(:,k);
    Q = Qres(:,:,k);
    sqrtQ = Q^(1/2);
    
    qNew = f(q, 0);
    
    x1 = @(v, q, sqrtQ) q(1) + sqrtQ(1,:)*v;
    x2 = @(v, q, sqrtQ) q(2) + sqrtQ(2,:)*v;

    x1_v = @(v) x1(v, q, sqrtQ);
    x2_v = @(v) x2(v, q, sqrtQ);

    g = @(args)...
        [ x1_v(args(1:2)) + dt*(args(3)+3)*x1_v(args(1:2))*(1 - x2_v(args(1:2)));...
          x2_v(args(1:2)) + dt*(args(3)+3)*x2_v(args(1:2))*(x1_v(args(1:2)) - 1)];

    cost1 = @(args) (-(g(args)-qNew)'*(g(args)-qNew));
    [opt1_1, fval1_1] = fmincon(cost1, [1;1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    [opt1_2, fval1_2] = fmincon(cost1, [1;-1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    [opt1_3, fval1_3] = fmincon(cost1, [-1;1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    [opt1_4, fval1_4] = fmincon(cost1, [-1;-1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    fval1s = [fval1_1, fval1_2, fval1_3, fval1_4];
    [~, idx] = min(fval1s);
    switch idx
        case 1
            opt1 = opt1_1;
        case 2
            opt1 = opt1_2;
        case 3
            opt1 = opt1_3;
        case 4
            opt1 = opt1_4;
    end
    p1 = g(opt1);
    r1 = norm(p1 - qNew);
    e1 = 1/r1 * (p1 - qNew);

    cost2 = @(args) (-norm( (g(args)-qNew) - e1*(e1'*(g(args)-qNew)) ));
    [opt2_1, fval2_1] = fmincon(cost2, [1;1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    [opt2_2, fval2_2] = fmincon(cost2, [1;-1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    [opt2_3, fval2_3] = fmincon(cost2, [-1;1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    [opt2_4, fval2_4] = fmincon(cost2, [-1;-1;0],...
        [], [],...% A, b
        [], [],...% Aeq, beq
        [-1; -1; -wMax], [1; 1; wMax],...% lb, ub
        @nonlincon, options);
    fval2s = [fval2_1, fval2_2, fval2_3, fval2_4];
    [~, idx] = min(fval2s);
    switch idx
        case 1
            opt2 = opt2_1;
        case 2
            opt2 = opt2_2;
        case 3
            opt2 = opt2_3;
        case 4
            opt2 = opt2_4;
    end
    p2 = g(opt2);
    e2 = (p2-qNew) - (e1*((p2-qNew)'*e1));
    e2 = e2 / norm(e2);
    r2 = norm(p2 - qNew);
    
    QNew = [e1, e2]*diag([r1,r2].^2)*[e1, e2]';

    qres(:,k+1) = qNew;
    Qres(:,:,k+1) = QNew;
end

function [c, ceq] = nonlincon(args)
c = norm(args(1:2)) - 1;
ceq = [];

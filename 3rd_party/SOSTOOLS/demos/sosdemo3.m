% SOSDEMO3 --- Bound on Global Extremum
% Section 3.3 of SOSTOOLS User's Manual
% 

clear; echo on;
syms x1 x2 gam;
vartable = [x1, x2];
% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);

% =============================================
% Declare decision variable gam too
prog = sosdecvar(prog,[gam]);

% =============================================
% Next, define SOSP constraints

% Constraint : r(x)*(f(x) - gam) >= 0
% f(x) is the Goldstein-Price function
% f1 = x1+x2+1;
% f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2;
% f3 = 2*x1-3*x2;
% f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2;

% f = (1+f1^2*f2)*(30+f3^2*f4);
f = 4*x1^2 - (21/10)*x1^4 + (1/3)*x1^6 + x1*x2 - 4*x2^2 + 4*x2^4;

prog = sosineq(prog,(f-gam));

% =============================================
% Set objective : maximize gam
prog = sossetobj(prog,-gam);

% =============================================
% And call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
% =============================================
% Finally, get solution
SOLgamma = sosgetsol(prog,gam)
echo off


%%

[x1, x2] = meshgrid(linspace(-1,1,51), linspace(-1,1,51));

fval = 4*x1.^2 - (21/10)*x1.^4 + (1/3)*x1.^6 + x1.*x2 - 4*x2.^2 + 4*x2.^4;

figure(321)
cla; hold on; grid on;
surf(x1,x2,fval)

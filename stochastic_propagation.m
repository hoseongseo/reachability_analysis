function P = stochastic_propagation(q, P_i, wMax, t)
% q: initial nominal state
% P_i: initial variance of each mixture
% wMax: maximum disturbance
% t: vector of time

dt = t(2) - t(1);
mu = q;

f = @(x,w) [(w+3)*x(1)*(1-x(2)); (w+3)*x(2)*(x(1)-1)]; % dynamics
g = @(x) [x(1)*(1-x(2)); x(2)*(x(1)-1)]; % f = 3*g + w*g
dgdx = @(x) [(1-x(2)), -x(1); x(2), (x(1)-1)];

wQ = wMax*(1/3);

% mean of mixtures at the initial
mu_1 = mu + 0.025*[1;0];
mu_2 = mu + 0.025*[0;1];
mu_3 = mu + 0.025*[-1;0];
mu_4 = mu + 0.025*[0;-1];

% trajectory of mean of mixtures
muTraj = zeros(2,4,length(t));
muTraj(:,:,1) = [mu_1, mu_2, mu_3, mu_4];
for k = 1:size(muTraj,2)
    for i = 1:length(t)-1
        x_ = muTraj(:,k,i);
        dx_ = f(x_,0);
        muTraj(:,k,i+1) = x_ + dt*dx_;
    end
end

wTraj = zeros(4,length(t));
wTraj(:,1) = 0.25*ones(4,1);

% uncertainty propagation
P = zeros(2,2,length(t));
P(:,:,1) = 0.25^2*4*P_i; % initial variance of sum of mixtures
P_1 = zeros(2,2,length(t));
P_1(:,:,1) = P_i;
P_2 = zeros(2,2,length(t));
P_2(:,:,1) = P_i;
P_3 = zeros(2,2,length(t));
P_3(:,:,1) = P_i;
P_4 = zeros(2,2,length(t));
P_4(:,:,1) = P_i;
for i = 1:length(t)-1
    P_1_ = P_1(:,:,i);
    mu_1_ = muTraj(:,1,i);
    A_1_ = 3*dgdx(mu_1_);
    g_1_ = g(mu_1_);
    dP_1_ = A_1_*P_1_ + P_1_*A_1_' + g_1_*(wQ^2)*g_1_';
    P_1(:,:,i+1) = P_1_ + dt*dP_1_;

    P_2_ = P_2(:,:,i);
    mu_2_ = muTraj(:,2,i);
    A_2_ = 3*dgdx(mu_2_);
    g_2_ = g(mu_2_);
    dP_2_ = A_2_*P_2_ + P_2_*A_2_' + g_2_*(wQ^2)*g_2_';
    P_2(:,:,i+1) = P_2_ + dt*dP_2_;

    P_3_ = P_3(:,:,i);
    mu_3_ = muTraj(:,3,i);
    A_3_ = 3*dgdx(mu_3_);
    g_3_ = g(mu_3_);
    dP_3_ = A_3_*P_3_ + P_3_*A_3_' + g_3_*(wQ^2)*g_3_';
    P_3(:,:,i+1) = P_3_ + dt*dP_3_;

    P_4_ = P_4(:,:,i);
    mu_4_ = muTraj(:,4,i);
    A_4_ = 3*dgdx(mu_4_);
    g_4_ = g(mu_4_);
    dP_4_ = A_4_*P_4_ + P_4_*A_4_' + g_4_*(wQ^2)*g_4_';
    P_4(:,:,i+1) = P_4_ + dt*dP_4_;
    
    mu_cur = muTraj(:,:,i);
    P_cur = [[P_1(1,1,i); P_1(1,2,i); P_1(2,2,i)],...
        [P_2(1,1,i); P_2(1,2,i); P_2(2,2,i)],...
        [P_3(1,1,i); P_3(1,2,i); P_3(2,2,i)],...
        [P_4(1,1,i); P_4(1,2,i); P_4(2,2,i)]];
    mu_next = muTraj(:,:,i+1);
    P_next = [[P_1(1,1,i+1); P_1(1,2,i+1); P_1(2,2,i+1)],...
        [P_2(1,1,i+1); P_2(1,2,i+1); P_2(2,2,i+1)],...
        [P_3(1,1,i+1); P_3(1,2,i+1); P_3(2,2,i+1)],...
        [P_4(1,1,i+1); P_4(1,2,i+1); P_4(2,2,i+1)]];
    
    [H1, H2, ~, ~] = coeff_integration(mu_cur, mu_next, P_cur, P_next, wQ,...
        [-0.1, 0.1], [201,201], dt);
    f2 = H2*wTraj(:,k);
	wOpt = quadprog(H1, f2, -eye(4), zeros(4,1), ones(1,4), 1, [], [], [],...
        optimset('display', 'none'));
    wTraj(:,i+1) = wOpt;
    
    P(:,:,i+1) = wTraj(1,i+1)^2*P_1(:,:,i+1)+...
        wTraj(2,i+1)^2*P_2(:,:,i+1)+...
        wTraj(3,i+1)^2*P_3(:,:,i+1)+...
        wTraj(4,i+1)^2*P_4(:,:,i+1);
end

function [H1, H2, x1gr, x2gr] =...
    coeff_integration(mu_cur, mu_next, P_cur, P_next, wMax, x_bin, n_bin, dt)

N = size(mu_cur, 2); % number of mixtures
mu_center = 0.5*(mean(mu_cur,2) + mean(mu_next,2));

[x1gr, x2gr] = meshgrid( mu_center(1) + linspace(x_bin(1), x_bin(2), n_bin(1)),...
    mu_center(2) + linspace(x_bin(1), x_bin(2), n_bin(2)) );
x_test = [x1gr(:), x2gr(:)]';
dx1 = x1gr(1,2)-x1gr(1,1);
dx2 = x2gr(2,1)-x2gr(1,1);

vals1 = zeros(prod(n_bin),N);
vals2 = zeros(prod(n_bin),N);
for k = 1:N
    vals1(:,k) = gaussian_2d(x_test, mu_next(:,k), P_next(:,k));
    vals2(:,k) = gaussian_2d(x_test, mu_cur(:,k), P_cur(:,k)) + dt*...
        lotka_volterra_coeff(x_test, wMax, mu_cur(:,k), P_cur(:,k));
end

H1 = zeros(N,N);
H2 = zeros(N,N);
for i = 1:N
    for j = 1:N
        H1(i,j) = sum(vals1(:,i).*vals1(:,j))*dx1*dx2;
        H2(i,j) = sum(vals1(:,i).*vals2(:,j))*dx1*dx2;
    end
end

function p = gaussian_2d(in1,in2,in3)
%GAUSSIAN_2D
%    P = GAUSSIAN_2D(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    26-Jul-2021 17:06:27

P1 = in3(1,:);
P2 = in3(2,:);
P3 = in3(3,:);
mu1 = in2(1,:);
mu2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
t2 = P2.^2;
t3 = pi.^2;
t4 = P1.*P3;
t5 = -t4;
t6 = t2+t5;
t7 = 1.0./t6;
p = exp(t7.*(mu2./2.0-x2./2.0).*(P1.*mu2-P2.*mu1-P1.*x2+P2.*x1)-t7.*(mu1./2.0-x1./2.0).*(P2.*mu2-P3.*mu1-P2.*x2+P3.*x1)).*1.0./sqrt(t2.*t3.*-4.0+t3.*t4.*4.0);

function cff = lotka_volterra_coeff(in1,wMax,in3,in4)
%LOTKA_VOLTERRA_COEFF
%    CFF = LOTKA_VOLTERRA_COEFF(IN1,WMAX,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    26-Jul-2021 17:02:51

P1 = in4(1,:);
P2 = in4(2,:);
P3 = in4(3,:);
mu1 = in3(1,:);
mu2 = in3(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
t2 = conj(wMax);
t3 = P1.*mu2;
t4 = P2.*mu1;
t5 = P2.*mu2;
t6 = P3.*mu1;
t7 = P1.*x2;
t8 = P2.*x1;
t9 = P2.*x2;
t10 = P3.*x1;
t11 = P2.^2;
t12 = x1.^2;
t13 = x2.^2;
t14 = P1.*P3;
t15 = x1.*x2.*2.0;
t16 = 1.0./pi;
t17 = x1-1.0;
t18 = x2-1.0;
t24 = mu1./2.0;
t25 = mu2./2.0;
t26 = x1./2.0;
t27 = x2./2.0;
t19 = -t4;
t20 = -t6;
t21 = -t7;
t22 = -t9;
t23 = -t15;
t28 = -t14;
t29 = -t11;
t30 = t17.^2;
t31 = t18.^2;
t32 = -t26;
t33 = -t27;
t34 = t11+t28;
t35 = t14+t29;
t36 = t24+t32;
t37 = t25+t33;
t40 = t3+t8+t19+t21;
t41 = t5+t10+t20+t22;
t38 = 1.0./t34;
t39 = 1.0./sqrt(t35);
t42 = t37.*t38.*t40;
t43 = t36.*t38.*t41;
t44 = -t43;
t45 = t42+t44;
t46 = exp(t45);
cff = (t16.*t39.*t46.*(wMax+x1.*3.0-x2.*3.0+(t12.*wMax)./2.0+(t13.*wMax)./2.0-(wMax.*x1)./2.0-(wMax.*x2)./2.0-wMax.*x1.*x2.*2.0))./2.0-(t16.*t39.*t46.*(t2.*t38.*t40.*x2.*(t12.*-2.0+t15+t17+x1-x2)-t2.*t38.*t41.*x1.*(t13.*-2.0+t15+t18-x1+x2)))./2.0-(t16.*t39.*t46.*(t7.*t27.*t30.*t38.*wMax+t10.*t26.*t31.*t38.*wMax+t8.*t17.*t18.*t38.*wMax.*x2))./2.0-(t16.*t39.*t46.*(t38.*t40.*(t17.*x2.*3.0+t2.*t27.*t30-(t2.*t18.*x1.*x2)./2.0)+t38.*t41.*(t18.*x1.*3.0-(t2.*t31.*x1)./2.0+t2.*t17.*t26.*x2)))./2.0+(t16.^2.*t38.^3.*wMax.*exp(t42.*2.0-t43.*2.0).*(t8.*t13+t3.*x2+t5.*x1+t10.*x1+t19.*x2+t20.*x1+t21.*x2-t3.*x1.*x2+t4.*x1.*x2-t5.*x1.*x2+t6.*x1.*x2+t7.*x1.*x2-t8.*x1.*x2-t10.*x1.*x2).^2)./8.0-(t16.*t39.*t46.*wMax.*(t12+t13-x1.*x2.*4.0+1.0))./2.0;

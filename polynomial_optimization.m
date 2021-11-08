clear; close all; clc;

addpath(genpath('3rd_party/helperOC-master')) % HJB equation solver
addpath(genpath('3rd_party/ToolboxLS'))
addpath(genpath('3rd_party/SOSTOOLS')) % SOS programming solver
addpath(genpath('3rd_party/SeDuMi_1_3')) % SDP solver (required for SOSTOOLS)

N = [4,4];
I = Polynomial.integer_grid(N);
c = randn( prod(N+1), 1 );
p = Polynomial(c, I );

B = inv(Polynomial.bernstein_transform_matrix(N));
h = prod((1./(I+1)).*( ones(size(I)).^(I+1) - (-ones(size(I))).^(I+1) ),2);

N_hat = [2,2];
I_hat = Polynomial.integer_grid(N_hat);
h_hat = prod((1./(I_hat+1)).*( ones(size(I_hat)).^(I_hat+1) - (-ones(size(I_hat))).^(I_hat+1) ),2);
E_hat = Polynomial.expand_matrix(N, I_hat);

%%
% without subdivision
cvx_begin
    variable c_hat1( prod(N_hat+1), 1 )
    
    cost = -h_hat'*c_hat1;
    minimize cost
%     minimize max( B*(c-E_hat*c_hat1) )

    subject to
    B * (c - E_hat * c_hat1) >= 0
%     B * (c - E_hat * c_hat1) <= 4.6
cvx_end
p_hat1 = Polynomial(c_hat1, I_hat);
disp(['cost = ', num2str(-h_hat'*c_hat1)])

%%
% with subdivision (sparse)
T11 = Polynomial.domain_transform_matrix( N, 0.5*ones(2,1), [-0.5; 0.5]);
T12 = Polynomial.domain_transform_matrix( N, 0.5*ones(2,1), [0.5; 0.5]);
T21 = Polynomial.domain_transform_matrix( N, 0.5*ones(2,1), [-0.5; -0.5]);
T22 = Polynomial.domain_transform_matrix( N, 0.5*ones(2,1), [0.5; -0.5]);
cvx_begin
    variable c_hat2( prod(N_hat+1), 1 )
    
    cost = -h_hat'*c_hat2;
    minimize cost
%     minimize max( B*(c-E_hat*c_hat2) )
    
    subject to
    B * T11 * (c - E_hat * c_hat2) >= 0
    B * T12 * (c - E_hat * c_hat2) >= 0
    B * T21 * (c - E_hat * c_hat2) >= 0
    B * T22 * (c - E_hat * c_hat2) >= 0
cvx_end
p_hat2 = Polynomial(c_hat2, I_hat);

%%
% with subdivision (fine)
T11 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-3;3]);
T12 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-1;3]);
T13 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[1;3]);
T14 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[3;3]);
T21 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-3;1]);
T22 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-1;1]);
T23 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[1;1]);
T24 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[3;1]);
T31 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-3;-1]);
T32 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-1;-1]);
T33 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[1;-1]);
T34 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[3;-1]);
T41 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-3;-3]);
T42 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[-1;-3]);
T43 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[1;-3]);
T44 = Polynomial.domain_transform_matrix( N, 0.25*ones(2,1), [0.25; 0.25].*[3;-3]);
cvx_begin
    variable c_hat3( prod(N_hat+1), 1 )
    
    cost = -h_hat'*c_hat3;
    minimize cost
%     minimize max( B*(c-E_hat*c_hat3) )

    subject to
    B * T11 * (c - E_hat * c_hat3) >= 0
    B * T12 * (c - E_hat * c_hat3) >= 0
    B * T13 * (c - E_hat * c_hat3) >= 0
    B * T14 * (c - E_hat * c_hat3) >= 0
    
    B * T21 * (c - E_hat * c_hat3) >= 0
    B * T22 * (c - E_hat * c_hat3) >= 0
    B * T23 * (c - E_hat * c_hat3) >= 0
    B * T24 * (c - E_hat * c_hat3) >= 0
    
    B * T31 * (c - E_hat * c_hat3) >= 0
    B * T32 * (c - E_hat * c_hat3) >= 0
    B * T33 * (c - E_hat * c_hat3) >= 0
    B * T34 * (c - E_hat * c_hat3) >= 0
    
    B * T41 * (c - E_hat * c_hat3) >= 0
    B * T42 * (c - E_hat * c_hat3) >= 0
    B * T43 * (c - E_hat * c_hat3) >= 0
    B * T44 * (c - E_hat * c_hat3) >= 0
cvx_end
p_hat3 = Polynomial(c_hat3, I_hat);

%%
% with subdivision (best)
% T11 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[-7;7]);
% T12 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[-5;7]);
% T13 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[-3;7]);
% T14 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[-1;7]);
% T15 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[1;7]);
% T16 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[3;7]);
% T17 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[5;7]);
% T18 = Polynomial.domain_transform_matrix( N, 0.125*ones(2,1), [0.125; 0.125].*[7;7]);

Ndiv = [5; 5];
argin = cell(1,length(Ndiv));
for j = 1:length(Ndiv)
    argin{j} = 1:Ndiv(j);
end
argout = cell(1,length(Ndiv));
[argout{:}] = ndgrid(argin{:});
Ndiv_gr = zeros(prod(Ndiv),length(Ndiv));
for j = 1:length(Ndiv)
    Ndiv_gr(:,j) = argout{j}(:);
end

cvx_begin
    variable c_hat4( prod(N_hat+1), 1 )
    
    cost = -h_hat'*c_hat4;
    minimize cost
%     minimize max( B*(c-E_hat*c_hat3) )

    subject to
    for k = 1:prod(Ndiv)
        Ndiv_gr_ = Ndiv_gr(k,:)'; 
        T_ = Polynomial.domain_transform_matrix( N, 1./Ndiv, (1./Ndiv).*( 2*(Ndiv_gr_-1) - (Ndiv-1) ) );
        B * T_ * (c - E_hat * c_hat4) >= 0
    end

cvx_end
p_hat4 = Polynomial(c_hat4, I_hat);

const = zeros(size(B,1), prod(Ndiv));
for k = 1:prod(Ndiv)
    Ndiv_gr_ = Ndiv_gr(k,:)';
    T_ = Polynomial.domain_transform_matrix( N, 1./Ndiv, (1./Ndiv).*( 2*(Ndiv_gr_-1) - (Ndiv-1) ) );
    
    const(:,k) = B * T_ * (c - E_hat * c_hat4);
end
%%
xgr = createGrid([-1,-1], [1,1], [51,51]);

figure(1)
cla; hold on; grid on; axis tight;

surf(xgr.xs{1}, xgr.xs{2}, p.eval( xgr ),...
    'linestyle', 'none', 'facealpha', 0.9)

% surf(xgr.xs{1}, xgr.xs{2}, p_hat1.eval( xgr ),...
%     'linestyle', '-', 'edgecolor', 'r',...
%     'facealpha', 0.0, 'facecolor', 'r')
% 
% surf(xgr.xs{1}, xgr.xs{2}, p_hat2.eval( xgr ),...
%     'linestyle', '-', 'edgecolor', 'b',...
%     'facealpha', 0.0, 'facecolor', 'b')
% 
% surf(xgr.xs{1}, xgr.xs{2}, p_hat3.eval( xgr ),...
%     'linestyle', '-', 'edgecolor', 'g',...
%     'facealpha', 0.0, 'facecolor', 'g')

surf(xgr.xs{1}, xgr.xs{2}, p_hat4.eval( xgr ),...
    'linestyle', '-', 'edgecolor', 'k',...
    'facealpha', 0.0, 'facecolor', 'k')
 
% surf(xgr.xs{1}, xgr.xs{2}, p.eval(xgr) - p_hat1.eval( xgr ),...
%     'linestyle', '-', 'edgecolor', 'r',...
%     'facealpha', 0.0, 'facecolor', 'r')
% 
% surf(xgr.xs{1}, xgr.xs{2}, p.eval(xgr) - p_hat2.eval( xgr ),...
%     'linestyle', '-', 'edgecolor', 'b',...
%     'facealpha', 0.0, 'facecolor', 'b')
% 
% surf(xgr.xs{1}, xgr.xs{2}, p.eval(xgr) - p_hat3.eval( xgr ),...
%     'linestyle', '-', 'edgecolor', 'g',...
%     'facealpha', 0.0, 'facecolor', 'g')
% 
% surf(xgr.xs{1}, xgr.xs{2}, p.eval(xgr) - p_hat4.eval( xgr ),...
%     'linestyle', '-', 'edgecolor', 'k',...
%     'facealpha', 0.0, 'facecolor', 'k')

xlabel('$x_1$')
ylabel('$x_2$')
% xlim([-1,0])
% ylim([0,1])
view([58,19])

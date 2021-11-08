clear; clc;

n = 3;

Hvec = sym('H', [n*(n+1)/2,1], 'real');
b = sym('b', [n,1], 'real');
r = sym('r', 'real');

D = Utils.duplication_matrix(n);

H = reshape( D*Hvec, [n,n] );
S = inv(H);
Svec = pinv(D)*S(:);

Q = (b'*inv(H)*b - r)*inv(H);
V = log(det(Q));
dVdc = jacobian(V, [r;b;Hvec]);


invQ = (1/(b'*inv(H)*b-r))*H;
dVdQ = invQ(:)'*D;

dQdxi_prime = [-Svec, 2*Svec*(b'*S), Svec*kron(b',b')*D + (b'*S*b - r)*eye(n*(n+1)/2)];
dxi_prime_dxi = blkdiag( eye(n+1), -pinv(D)*kron(S,S)*D );
dVdc_tmp = dVdQ * dQdxi_prime * dxi_prime_dxi;

for i = 1:length(dVdc)
    disp(['i = ', num2str(i), ' / ', num2str(length(dVdc))])
    simplify(dVdc_tmp(i) - dVdc(i))
end

%%


Q2 = (b'*S*b - r)*S;
V2 = log(det(Q2));
dV2dc = jacobian(V, [r;b;Svec]);

%%
vecQ = [Q(:,1); Q(2:end,2); Q(3,3)];
dQdc = jacobian(vecQ, [r;b;Hvec]);

vecQ2 = [Q2(:,1); Q2(2:end,2); Q2(3,3)];
dQ2dc = jacobian(vecQ2, [r;b;Svec]);

%%
dVdc_tmp = dVdQ'*dQdc;

%%
xi = [r; b; Hvec];
invH = inv(H);
Svec = pinv(D)*invH(:);
xi_prime = [r; b; Svec];
J = jacobian(xi_prime, xi);

M = zeros(n*n,n*n);
cnt = 1;
for i = 1:n
    e_i = zeros(n,1);
    e_i(i) = 1.0;
    for j = 1:n
        e_j = zeros(n,1);
        e_j(j) = 1.0;
        tmp_ = e_i * e_j';
        M(:,cnt) = tmp_(:);
        cnt = cnt + 1;
    end
end

ttt = simplify(-pinv(D) * kron(invH,invH) * D)
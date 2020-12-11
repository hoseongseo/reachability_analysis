function S = dyn_coeff(fp, n, m, N)
% fp : polynomial dynamics
% n : # of state
% m : # of disturbance
if nargin < 4
    order_all = reshape([fp.order], [numel([fp.order])/(n+m), (n+m)]);
    N = max(order_all,[],1);
end
Nx = N(1:n);
Nw = N(n+(1:m));

S = zeros( prod(Nx+1), prod(Nw+1), n );
for i = 1:n
    fpw = fp(i).collect(n+(1:m));
    for j = 1:prod(Nw+1)
        S(:,j,i) = Polynomial.expand_matrix(Nx, fpw.coeff(j).order)*fpw.coeff(j).coeff;
    end
end
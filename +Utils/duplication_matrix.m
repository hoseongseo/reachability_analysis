function D = duplication_matrix(n)
% D*vech(A) = vec(A)
D = zeros(n*n, 0.5*n*(n+1));
for j = 1:n
    for i = j:n
        u_ = zeros(0.5*n*(n+1),1);
        u_( (j-1)*n + i - 0.5*j*(j-1) ) = 1;
        T_ = zeros(n,n);
        T_(i,j) = 1;
        T_(j,i) = 1;
        D = D + T_(:)*u_';
    end
end
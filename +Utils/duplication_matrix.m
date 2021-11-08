function D = duplication_matrix(n, order)
% D*vech(A) = vec(A)
D = zeros(n*n, 0.5*n*(n+1));
if nargin < 2
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
else
    % compute idx
    idx = zeros(size(order,1), 2);
    for i = 1:size(order,1)
        order_ = order(i,:);
        cnt = 1;
        for k = 1:size(order,2)
            if order_(k) > 1
                idx(i,1) = k;
                idx(i,2) = k;
                break;
            elseif order_(k) > 0
                idx(i,cnt) = k;
                cnt = cnt + 1;
            end
            if cnt > 2
                break;
            end
        end
    end
    
    for k = 1:size(idx,1)
        i = idx(k,1);
        j = idx(k,2);
        u_ = zeros(0.5*n*(n+1),1);
        u_( (i-1)*n + j - 0.5*i*(i-1) ) = 1;
        T_ = zeros(n,n);
        T_(i,j) = 1;
        T_(j,i) = 1;
        D = D + T_(:)*u_';
    end
end
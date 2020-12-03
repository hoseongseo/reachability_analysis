function D = derivative_matrix(N)
n = length(N);
I = integer_grid(zeros(size(N)), N);

D = [];
for i = 1:n
    e = zeros(1,n);
    e(i) = 1;
    
    I_ = I - repmat(e, [prod(N+1),1]);
    idx_ = I_(:,i) >= 0;
    I_ = I_(idx_,:);
    
    A_ = zeros(sum(idx_), length(idx_));
    A_(:,idx_) = diag(I(idx_,i));
    D_ = Polynomial.expand_matrix(N, I_) * A_;

    D = cat(3, D, D_);
end
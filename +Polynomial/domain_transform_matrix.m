function T = domain_transform_matrix(n, a, b)
I = Polynomial.integer_grid(zeros(size(n)), n);
T = zeros( prod(n+1), prod(n+1) );
if class(a) == "sym" || class(b) == "sym"
    T = sym(T);
end
    
for i = 1:prod(n+1)
    i_ = I(i,:);
    K_ = Polynomial.integer_grid(i_, n);
    id_K_ = Polynomial.getId( K_, n );
    for k = 1:size(K_,1)
        k_ = K_(k,:);
        val_ = 1;
        for j = 1:length(n)
            val_ = val_ * nchoosek(k_(j),i_(j)) * a(j)^i_(j) * b(j)^(k_(j)-i_(j));
        end
        T(i, id_K_(k)) = val_;
    end
end
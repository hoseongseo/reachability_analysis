function [S, Bp] = bernstein_transform_matrix(n)
I = integer_grid(zeros(size(n)), n);
vp(1:length(n)) = Polynomial;
for i = 1:length(n)
    e_ = zeros(1,length(n));
    e_(i) = 1.0;
    vp(i) = Polynomial(1, e_);
end

Bp(1:prod(n+1)) = Polynomial(1, zeros(1,length(n)));
for i = 1:prod(n+1)
    for j = 1:length(n)
        Bp(i) = Bp(i) * nchoosek(n(j), I(i,j)) * ...
            ( (vp(j)+1)*0.5 ).^I(i,j) *...
            ( (1-vp(j))*0.5 ).^(n(j)-I(i,j));
    end
end
S = [Bp.coeff];
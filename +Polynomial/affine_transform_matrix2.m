function [T, full_order, x_prods] = affine_transform_matrix2(I, A, b)
% [x_1, x_1^2, x_1^3, ...]
% [x_2, x_2^2, x_2^3, ...]
% [ ...
% [x_n, x_n^2, x_n^3, ...]
n = max(I, [], 1);
v_order = [zeros(size(n)); eye(length(n))];
x_powers(1:length(n), 1:max(n,[],2)+1) = Polynomial(1, zeros(size(n)));
for j = 1:length(n)
    for i = 1:n(j)
        x_powers(j,i+1) = x_powers(j,i) .* Polynomial([b(j); A(j,:)'], v_order);
    end
end

%I = Polynomial.integer_grid(zeros(size(n)), n);
x_prods = x_powers(1, I(:,1)+1);
for j = 2:length(n)
    x_prods = x_prods .* x_powers(j, I(:,j)+1);
end

% max_order = max(x_prods(end).order,[],1);

full_order = x_prods(end).order;
full_order_sum = sum(full_order,2);
% T = zeros(size(full_order,1), prod(n+1));
T = zeros(size(full_order,1), size(I,1));
if class(b) == "sym"
    T = sym(T);
end
for k = 1:size(I,1)
    cur_sum = x_prods(k).order(end,end);
    cur_idx = full_order_sum <= cur_sum;
    C = zeros(size(full_order,1), length(x_prods(k).coeff));
    C( cur_idx, : ) = eye(sum(cur_idx));
    T(:,k) = C * x_prods(k).coeff;
end


% max_order = sum(n,2) * ones(1,size(n,2)); % possible maximum order
% I = integer_grid(zeros(size(n)), max_order);

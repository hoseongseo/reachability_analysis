function C = expand_matrix(max_order, order)
id = Polynomial.getId(order, max_order);
C = zeros(prod(max_order+1), length(id));
C(id,:) = eye(length(id));
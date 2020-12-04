function P = product_matrix(N1, N2)

I1 = Polynomial.integer_grid(zeros(size(N1)), N1);
I2 = Polynomial.integer_grid(zeros(size(N2)), N2);
N = N1 + N2;
id1 = Polynomial.getId(I1, N);
id2 = Polynomial.getId(I2, N);

id_all = Utils.minksum(id1, id2)-1;
id = unique(id_all);
P = zeros(length(id), length(id_all));
for i = 1:length(id)
    P(i,:) = (id_all == id(i));
end

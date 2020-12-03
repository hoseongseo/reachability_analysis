function Hvec = shape2vec(H)
n = size(H,1);
L = n*(n+1)/2;
Hvec = zeros(L,1);
for i = 1:n
    Hvec( (i-1)*(i)/2 + (1:i) ) = H(1:i,i);
end
end
function i = integer_grid(n1, n2)

m = length(n1);
argin = cell(1,m);
argout = cell(1,m);

if nargin > 1
    for j = 1:m
        argin{j} = n1(j):n2(j);
    end
else
    for j = 1:m
        argin{j} = 0:n1(j);
    end
end
[argout{:}] = ndgrid( argin{:} );

if nargin > 1
    i = zeros(prod(n2-n1+1),m);
else
    i = zeros(prod(n1+1),m);
end

for j = 1:m
    i(:,j) = argout{j}(:);
end
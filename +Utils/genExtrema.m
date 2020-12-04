function pts = genExtrema(bnd)
% generate minmax points
N = length(bnd);
d = (0:(2^N - 1))';
result = [];
for i = 0:(N - 1)
    result = [result, round(rem(d * (2^i), (2^(N))) / (2^(N)))];
end
pm = (result - (1 - result))';
pts = pm.*repmat(bnd, [1,size(pm,2)]);

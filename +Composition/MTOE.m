function [Qsum, alpha] = MTOE(Q)

N = size(Q,3);
tr = zeros(1,N);
for k = 1:N
    tr_ = trace(Q(:,:,k));
    tr(k) = sqrt(tr_);
end
Qsum = zeros(size(Q(:,:,1)));
for k = 1:N
    Qsum = Qsum + Q(:,:,k) / tr(k);
end
Qsum = Qsum * sum(tr);
alpha = tr / sum(tr);
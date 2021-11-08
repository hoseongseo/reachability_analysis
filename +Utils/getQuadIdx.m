function idx = getQuadIdx(order)
idx = zeros(size(order,1), 2);
for i = 1:size(order,1)
    order_ = order(i,:);
    cnt = 1;
    for k = 1:size(order,2)
        if order_(k) > 1
            idx(i,1) = k;
            idx(i,2) = k;
            break;
        elseif order_(k) > 0
            idx(i,cnt) = k;
            cnt = cnt + 1;
        end
        if cnt > 2
            break;
        end
    end
end
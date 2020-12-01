function pts = getLevelSet(gr, V, val)
% gr : grid structure
% V : value function defined on the grid
% val : target level

if gr.dim == 1
    pts = getLevelSet1d(gr, V, val);
elseif gr.dim == 2
    pts = getLevelSet2d(gr, V, val);
elseif gr.dim == 3
    pts = getLevelSet3d(gr, V, val);
else
    error('Dimension larger than 3 are not supported!!')
end
end

function pts = getLevelSet1d(gr, V, val)
v = gr.vs{1};
zci = @(v) find((v(:)-val).*circshift(v(:)-val, [-1 0]) <= 0);
bdry_idx = zci( V );
pts = zeros(1,length(bdry_idx));
for i = 1:length(bdry_idx)
    tmp_idx_ = bdry_idx(i);
    v1_ = v(tmp_idx_);
    V1_ = V(tmp_idx_);
    v2_ = v(tmp_idx_-sign(V1_));
    V2_ = V(tmp_idx_-sign(V1_));
    v_ = v1_ + (-V1_)/(V2_-V1_)*(v2_-v1_);
    pts(i) = v_;
end
% this function sometimes results error
if length(pts) < 2
    if V(1) < 0
        pts = [v(1), pts];
    else
        pts = [pts, v(end)];
    end
end
end

function pts = getLevelSet2d(gr, V, val)
v1 = gr.vs{1};
v2 = gr.vs{2};
slice = contourc(v1, v2, V', [val,val]);
pts = [];
if ~isempty(slice)
    C_ind = find(slice(1,:) == val);
    N_contour = slice(2, C_ind);
%     [~, max_ind] = max(N_contour);
%     slice = slice(:,C_ind(max_ind)+[1:N_contour(max_ind)]);
%     pts = slice;
    for i = 1:length(N_contour)
        if N_contour(i) > 1
            pts = [pts, slice(:,C_ind(i)+[1:N_contour(i)])];
        end
    end
end
end

function pts = getLevelSet3d(gr, V, val)
v1 = gr.vs{1};
v2 = gr.vs{2};
v3 = gr.vs{3};
pts = [];
for i = 1:length(v3) % along z axis
    slice = contourc(v1, v2, V(:,:,i)', [val,val]);
    if ~isempty(slice)
        C_ind = find(slice(1,:) == val);
        N_contour = slice(2, C_ind);
        for k = 1:length(N_contour)
            if N_contour(k) > 1
                slice2 = slice(:,C_ind(k)+[1:N_contour(k)]);
                pts = [pts, [slice2; v3(i)*ones(1,size(slice2,2))]];
            end
        end
%         [~, max_ind] = max(N_contour);
%         slice = slice(:,C_ind(max_ind)+[1:N_contour(max_ind)]);
%         pts = [pts, [slice; v3(i)*ones(1,size(slice,2))]];
    end
end
for i = 1:length(v1) % along x axis
    slice = contourc(v2, v3, squeeze(V(i,:,:))', [val,val]);
    if ~isempty(slice)
        C_ind = find(slice(1,:) == val);
        N_contour = slice(2, C_ind);
        for k = 1:length(N_contour)
            if N_contour(k) > 1
                slice2 = slice(:,C_ind(k)+[1:N_contour(k)]);
                pts = [pts, [v1(i)*ones(1,size(slice2,2)); slice2]];
            end
        end
%         [~, max_ind] = max(N_contour);
%         slice = slice(:,C_ind(max_ind)+[1:N_contour(max_ind)]);
%         pts = [pts, [v1(i)*ones(1,size(slice,2)); slice]];
    end
end
for i = 1:length(v2) % along y axis
    slice = contourc(v1, v3, squeeze(V(:,i,:))', [val,val]);
    if ~isempty(slice)
        C_ind = find(slice(1,:) == val);
        N_contour = slice(2, C_ind);
        for k = 1:length(N_contour)
            if N_contour(k) > 1
                slice2 = slice(:,C_ind(k)+[1:N_contour(k)]);
                pts = [pts, [slice2(1,:); v2(i)*ones(1,size(slice2,2)); slice2(2,:)]];
            end
        end
%         [~, max_ind] = max(N_contour);
%         slice = slice(:,C_ind(max_ind)+[1:N_contour(max_ind)]);
%         pts = [pts, [slice(1,:); v2(i)*ones(1,size(slice,2)); slice(2,:)]];
    end
end
end

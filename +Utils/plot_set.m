function hh = plot_set(bdry, t, n, color, alpha)

for i = 1:round(length(t)/n):length(t)
    if class(bdry) == "cell"
        pts = bdry{i};
    else
        pts = bdry(:,:,i);
    end
    
    pts = [pts(1,:); t(i)*ones(1,size(pts,2)); pts(2,:)]; 
    if i == 1
        pts_prev = pts;
    else
        pts_total = unique([pts_prev, pts]', 'rows')';

        hi = find(pts_total(2,:) == pts(2,1));
        K = convhulln(pts_total', {'Qt', 'Pp'});
        K_new = [];
        for k = 1:size(K,1)
            tf = ismember(K(k,:), hi);
            if sum(tf) > 0 && sum(tf) < 3
                K_new = [K_new; K(k,:)];
            end
        end
        
        hh = trisurf(K_new, pts_total(1,:), pts_total(2,:), pts_total(3,:),...
            'facecolor', color, 'facealpha', alpha, 'linestyle', 'none', 'linewidth', 0.01);
        
        pts_prev = pts;     
    end    
    disp(['Plotting i = ', num2str(i), ' / ', num2str(length(t))])
end
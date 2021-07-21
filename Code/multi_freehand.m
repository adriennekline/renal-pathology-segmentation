function multi_mask = multi_freehand(num_roi,image_size)
% multiple ROI selection
% Adrienne Kline
% University of Calgary
% Copyright (c) 2020

    for i = 1:num_roi
        h(i) = drawfreehand;
        h(i).FaceAlpha = 1;
        h(i).FaceSelectable = false;
        freehanding{i} = createMask(h(i),image_size);
    end
    multi_mask = zeros(size(image_size));
    for j = 1:size(freehanding,2)
        multi_mask = multi_mask + freehanding{j};
    end
    clear j i
end


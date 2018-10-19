function [x,flip] = ROIcrop(px, py, xlim)

if px(1) > px(2)
    tx = px(1);
    ty = py(1);
    px(1) = px(2);
    py(1) = py(2);
    px(2) = tx;
    py(2) = ty;
end

m = (py(2)-py(1))/(px(2)-px(1));

if m > 0
    flip = 0;
    x(1) = max(px(1)-40,1);
    x(3) = min(px(2)+15,xlim);
    x(2) = py(1);
    x(4) = py(2);
else
    flip = 1;
    x(1) = max(px(1)-15,1);
    x(3) = min(px(2)+40,xlim);
    x(2) = py(2);
    x(4) = py(1);
end

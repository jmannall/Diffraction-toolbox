function [k, validSource] = PointPlanePosition(source, planeNormal, d)
    
    k = dot(source, planeNormal) - d;
    if k < 0
        validSource = false;
    else
        validSource = true;
    end
end
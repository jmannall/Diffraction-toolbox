function [k, validSource] = PointPlanePosition(source, planeNormal, d)
    
    % k is the distance along the normal from the plane?
    k = dot(source, planeNormal) - d;   % Equation of a plane: ax + by + cz = d. where (a, b, z) are plane normal. If source on plane k = 0. If k > 0 source in front of plane etc
    if k <= 0
        validSource = false;
    else
        validSource = true;
    end
end
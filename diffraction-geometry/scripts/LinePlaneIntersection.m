function [intersection, valid] = LinePlaneIntersection(lineStart, lineEnd, normal, d)

    gradient = (lineStart - lineEnd);
    scale = dot(normal, gradient);
    k = dot(lineStart, normal) - d;
    intersection = lineStart - gradient .* k / scale;

    % Check intersection lies within line sector
    gradient = (lineStart - intersection);
    scaleI = dot(normal, gradient);
    if scaleI > scale && scaleI > 0
        valid = false;
    elseif scaleI < scale && scaleI < 0
        valid = false;
    else
        valid = true;
    end
end
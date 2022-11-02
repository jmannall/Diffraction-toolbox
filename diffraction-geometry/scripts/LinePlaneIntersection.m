function intersection = LinePlaneIntersection(lineStart, lineEnd, normal, d)

    gradient = (lineStart - lineEnd);
    scale = dot(normal, gradient);
    k = dot(lineStart, normal) - d;
    intersection = lineStart - gradient .* k / scale;
end
function [intersection, validPath] = CheckValidLinePlaneIntersection(lineStart, lineEnd, normal, d, planes, corners, numPlaneCorners, i)
    intersection = LinePlaneIntersection(lineStart, lineEnd, normal, d);

    validCorners = planes(i,:) > 0;
    planeCorners = [corners(planes(i,validCorners),:); corners(planes(i,1),:)];
    angleRotation = zeros(numPlaneCorners(i), 1);
    for j = 1:numPlaneCorners(i)
        vecOne = intersection - planeCorners(j,:);
        vecTwo = intersection - planeCorners(j + 1,:);
        dotProduct = dot(normal, cross(vecOne, vecTwo));
        angleRotation(j) = sign(dotProduct) * acosd(dot(vecOne, vecTwo) / (norm(vecOne) * norm(vecTwo)));
    end
    if round(sum(angleRotation), 5) ~= 360
        validPath = false;
    else
        validPath = true;
    end
end
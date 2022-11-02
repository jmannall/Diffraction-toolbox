function obstruction = CheckForObstruction(source, receiver, planes, normals, d, corners, numPlaneCorners, numPlanes, i, obstruction)
    
    j = 1;
    while ~obstruction && j <= numPlanes
            if j ~= i
                [kS, ~] = PointPlanePosition(source, normals(j,:), d(j));
                [kR, ~] = PointPlanePosition(receiver, normals(j,:), d(j));
                if kS * kR < 0
                    [~, obstruction] = CheckValidLinePlaneIntersection(receiver, source, normals(j,:), d(j), planes, corners, numPlaneCorners, j);
                end
            end
            j = j + 1;
    end
end
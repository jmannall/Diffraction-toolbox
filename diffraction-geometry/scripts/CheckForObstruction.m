function obstruction = CheckForObstruction(source, receiver, room, i, obstruction)

    if nargin < 5
        obstruction = false;
    end
    j = 1;
    while ~obstruction && j <= room.numPlanes
            if j ~= i
                [kS, ~] = PointPlanePosition(source, room.planeNormals(j,:), room.d(j));
                [kR, ~] = PointPlanePosition(receiver, room.planeNormals(j,:), room.d(j));
                if kS * kR < 0
                    [~, obstruction] = CheckValidLinePlaneIntersection(receiver, source, room, j);
                end
            end
            j = j + 1;
    end
end
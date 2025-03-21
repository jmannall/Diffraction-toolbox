function [intersection, validPath] = CheckValidLinePlaneIntersection(lineStart, lineEnd, room, i)

    [intersection, validPath] = LinePlaneIntersection(lineStart, lineEnd, room.planeNormals(i,:), room.d(i));

    if validPath
        validCorners = room.planeCorners(i,:) > 0;
        planeCorners = [room.corners(room.planeCorners(i,validCorners),:); room.corners(room.planeCorners(i,1),:)];
        angleRotation = zeros(room.numCorners(i), 1);
        for j = 1:room.numCorners(i)
            vecOne = intersection - planeCorners(j,:);
            vecTwo = intersection - planeCorners(j + 1,:);
            dotProduct = dot(room.planeNormals(i,:), cross(vecOne, vecTwo));
            angleRotation(j) = sign(dotProduct) * acosd(dot(vecOne, vecTwo) / (norm(vecOne) * norm(vecTwo)));
        end
        if round(sum(angleRotation), 5) ~= 360
            validPath = false;
        end
    end
end
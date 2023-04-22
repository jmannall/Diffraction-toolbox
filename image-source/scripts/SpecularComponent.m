function [specular, plot, sData] = SpecularComponent(room, plot)

disp('Reflection order: 1')
idx = 0;
vIdx = 0;
specular.valid = false;

tic
for i = 1:room.numPlanes
    disp(['Plane: ', num2str(i)])
    % Check if source behind the plane
    [k, validSource] = PointPlanePosition(room.source, room.planeNormals(i,:), room.d(i));

    if validSource
        % Create virtual source
        vS = room.source - 2 * room.planeNormals(i,:) * k;
        pathLength = PathLength(vS, room.receiver);
        if pathLength < room.maxPathLength
            idx = idx + 1;  % Source will be stored
            sData.vS(idx,:) = vS;
            sData.planes(idx,:) = i;
            sData.idxs(idx,:) = idx;
    
            % Find intersection point and check it is valid
            [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.receiver, vS, room, i);
    
            if validIntersection && room.receiverCanSeePlane(i)
                % Check if path is blocked
                obstruction = CheckForObstruction(intersection, room.receiver, room, i);
                obstruction = CheckForObstruction(room.source, intersection, room, i, obstruction);
                if ~obstruction
                    specular.valid = true;
                    vIdx = vIdx + 1;
                    specular.pathLength(vIdx) = pathLength;
                    specular.planes(vIdx) = i;
                    plot.specular{vIdx} = [room.source; intersection; room.receiver];
                end
            end
        end
    end
    sData.n = idx;
end
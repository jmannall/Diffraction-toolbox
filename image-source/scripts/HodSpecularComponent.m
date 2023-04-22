function [hodSpecular, plot, sData] = HodSpecularComponent(room, plot, data, refOrder)

    hodSpecular = cell(refOrder, 1);
    sData = cell(refOrder, 1);
    
    sData{1} = data;
    for j = 2:refOrder
        specular.valid = false;
        idx = 0;
        vIdx = 0;
    
        disp(['Reflection order: ', num2str(j)])
    
        tic
        for i = 1:room.numPlanes
            disp(['Plane: ', num2str(i)])
            for n = 1:data.n
                lastPlane = data.planes(n);
                if lastPlane ~= i && room.planeCanSeeplane(i, lastPlane)
                    % Check if source behind the plane
                    [k, validSource] = PointPlanePosition(data.vS(n,:), room.planeNormals(i,:), room.d(i));
                    
                    if validSource
                        % Create virtual source
                        vS = data.vS(n,:) - 2 * room.planeNormals(i,:) * k;
                        pathLength = norm(vS - room.receiver);
                        if pathLength < room.maxPathLength
                            idx = idx + 1;  % Source will be stored
                            sData{j}.vS(idx,:) = vS;
                            sData{j}.planes(idx,:) = [data.planes(n,:), i];
                            sData{j}.idxs(idx,:) = [data.idxs(n,:), idx];
                            
                            % Find intersection point and check it is valid
                            intersection = zeros(j, 3);
                            [intersection(j,:), validIntersection] = CheckValidLinePlaneIntersection(room.receiver, vS, room, i);
                            m = j - 1;
                            while validIntersection && m > 0
                                refPlane = data.planes(n,m);
                                startPos = sData{m}.vS(data.idxs(n,m),:);
                                [intersection(m,:), validIntersection] = CheckValidLinePlaneIntersection(intersection(m + 1,:), startPos, room, refPlane);
                                m = m - 1;
                            end
    
                            if validIntersection && room.receiverCanSeePlane(i)
                                % Check if path is blocked
                                obstruction = CheckForObstruction(intersection(j,:), room.receiver, room, i);
                                obstruction = CheckForObstruction(room.source, intersection(1,:), room, i, obstruction);
                                p = 1;
                                while ~obstruction && p < j
                                    obstruction = CheckForObstruction(intersection(p,:), intersection(p + 1,:), room, i, obstruction);
                                    p = p + 1;
                                end
                                if ~obstruction
                                    specular.valid = true;
                                    vIdx = vIdx + 1;
                                    specular.pathLength(vIdx) = pathLength;
                                    specular.planes(vIdx) = i;
                                    plot.hodSpecular{j}{vIdx} = [room.source; intersection; room.receiver];
                                end
                            end
                        end
                    end
                end
            end
        end
        toc
        hodSpecular{j} = specular;
        clear specular
        sData{j}.n = idx;
        data = sData{j};
    end
end
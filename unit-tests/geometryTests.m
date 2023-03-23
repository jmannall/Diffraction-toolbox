
%% Shoebox room
x = 2;
y = 3;
z = 2.6;

[room, source, receiver] = CreateShoeboxRoomGeometry(x, y, z);
PlotGeometry(room.corners, room.planeCorners, source, receiver)

%% L shaped room

x = [2 4];
y = [3 6];
z = 2.6;

[room, source, receiver] = CreateLShapedRoomGeometry(x, y, z);
PlotGeometry(room.corners, room.planeCorners, source, receiver)

%% Ring shaped room

close all
clear all

x = [8 2 4];
y = [12 3 6];
z = 2.6;

[room, source, receiver] = CreateRingShapedRoomGeometry(x, y, z);
%source = [6, 2, 1.6];

PlotGeometry(room.corners, room.planeCorners, source, receiver)

refOrder = 5;
diffOrder = 3;
combinedOrder = 5;
maxPathLength = 100;

room.numPlanes = size(room.planeCorners, 1);
room.numEdges = size(room.edgeCorners, 1);
room.planeNormals = zeros(room.numPlanes, 3);
[room.d, room.numCorners, kR] = deal(zeros(room.numPlanes, 1));
receiverCanSeePlane = true(room.numPlanes, 1);

%% Create geometry parameters

for i = 1:room.numPlanes
    validCorners = room.planeCorners(i,:) > 0; % Find non zero plane corners
    room.planeCornersCoords{i} = [room.corners(room.planeCorners(i,validCorners),:); room.corners(room.planeCorners(i,1),:)];
    room.numCorners(i) = sum(validCorners);     % num corners in the plane
    edge1 = room.corners(room.planeCorners(i,1),:) - room.corners(room.planeCorners(i,2),:);
    edge2 = room.corners(room.planeCorners(i,1),:) - room.corners(room.planeCorners(i,3),:);
    room.planeNormals(i,:) = normr(cross(edge1, edge2));    % Plane normals
    room.d(i) = dot(room.planeNormals(i,:), room.corners(room.planeCorners(i,1),:));     % Equation of a plane: ax + by + cz = d. where (a, b, z) are plane normal
    % Check if receiver lies behind the plane
    [kR(i), receiverCanSeePlane(i)] = PointPlanePosition(receiver, room.planeNormals(i,:), room.d(i));
end

% Source and receiver to edge visibility (assuming no partially obstructed
% edges)
[sourceCanSeeEdge, receiverCanSeeEdge] = deal(false(room.numEdges, 1));
for i = 1:room.numEdges
    room.edgeMidPoint(i,:) = (room.corners(room.edgeCorners(i,1),:) + room.corners(room.edgeCorners(i,2),:)) / 2;
    receiverCanSeeEdge(i) = ~CheckForObstruction(room.edgeMidPoint(i,:), receiver, room, 0);
    sourceCanSeeEdge(i) = ~CheckForObstruction(source, room.edgeMidPoint(i,:), room, 0);
end

planeCanSeeplane = true(room.numPlanes);
edgeCanSeePlane = false(room.numEdges, room.numPlanes);
for i = 1:room.numPlanes
    for j = 1:room.numPlanes
        dotProduct = dot(room.planeNormals(i,:), room.planeNormals(j,:));
        if dotProduct == 1  % Normals are parallel and in the same direction
            [planeCanSeeplane(i,j), planeCanSeeplane(j,i)] = deal(false);
        else
            validCorners = false(1, room.numCorners(j));
            for n = 1:room.numCorners(j)
                [~, validCorners(n)] = PointPlanePosition(room.corners(room.planeCorners(j,n),:), room.planeNormals(i,:), room.d(i));  % Check if each plane corner is behind other the plane
            end
            if sum(validCorners) == 0
                [planeCanSeeplane(i,j), planeCanSeeplane(j,i)] = deal(false);   % If all corners are behind the other plane. No reflection is every possible between them
            end
        end
    end
    for j = 1:room.numEdges
        [~, edgeCanSeePlane(j,i)] = PointPlanePosition(room.edgeMidPoint(j,:), room.planeNormals(i,:), room.d(i));  % Check if each plane corner is behind other the plane
    end
end

%% Direct sound

i = 0;
lineOfSight = ~CheckForObstruction(source, receiver, room, i);

if lineOfSight
    directPathLength = PathLength(receiver, source);
end

%% First order diffraction

inf = 1e5;
diffPath = zeros(inf, diffOrder);
diffIdx = 0;
edgePath = {};

tic
disp('Diffraction order: 1')
for i = 1:room.numEdges
    if receiverCanSeeEdge(i) && sourceCanSeeEdge(i)
        diffIdx = diffIdx + 1;
        validDiffPath(diffIdx) = true;
        diffPath(diffIdx,:) = [i, zeros(1, diffOrder - 1)];
        edgePath{diffIdx} = [source; room.edgeMidPoint(i,:); receiver];
    end
end
toc

% PlotGeometry(room.corners, room.planeCorners, source, receiver, lineOfSight, vSources, intersections, validPath, edgePath)

%% Higher order diffraction

edgeStore = 1:room.numEdges;
edges = edgeStore;
for i = 1:diffOrder
    edges = [edges, edgeStore];
end

for j = 2:diffOrder
    tic
    disp(['Diffraction order: ', num2str(j)])
    paths = unique(nchoosek(edges,j), "rows");
    numPaths = size(paths, 1);
    for i = 1:numPaths
        path = paths(i,:);
        if sourceCanSeeEdge(path(1)) && receiverCanSeeEdge(path(j))
            valid = true;
            n = 1;
            while valid && n < j
                valid = room.edgeCanSeeEdge(path(n), path(n + 1));
                n = n + 1;
            end
            if valid
                diffIdx = diffIdx + 1;
                validDiffPath(diffIdx) = true;
                diffPath(diffIdx,:) = [path, zeros(1, diffOrder - j)];
                edgePath{diffIdx} = [source; room.edgeMidPoint(path,:); receiver];
            end
        end
    end
    toc
end

% PlotGeometry(room.corners, room.planeCorners, source, receiver, lineOfSight, vSources, intersections, validPath, edgePath)

%% First order reflections

vSources = zeros(inf, 3);
validPath = false(inf, 1);
[pathIdx, specularPathLength, refOrderPaths] = deal(zeros(inf, 1));
[refPlanePath, vSourceIdxPath] = deal(zeros(inf, refOrder));

idx = 0;
count = 0;
tic
disp('Reflection order: 1')
for i = 1:room.numPlanes
    disp(['Plane: ', num2str(i)])
    % Check if source behind the plane
    [k, validSource] = PointPlanePosition(source, room.planeNormals(i,:), room.d(i));

    if validSource
        % Create virtual source
        vS = source - 2 * room.planeNormals(i,:) * k;
        pathLength = PathLength(vS, receiver);
        if pathLength > maxPathLength
            validSource = false;
        end
        if validSource
            idx = idx + 1;  % Source will be stored
            vSources(idx,:) = vS;
            refPlanePath(idx,:) = [i, zeros(1, refOrder - 1)];
            vSourceIdxPath(idx,:) = [idx, zeros(1, refOrder - 1)];
    
            % Find intersection point and check it is valid
            [intersection, validIntersection] = CheckValidLinePlaneIntersection(receiver, vSources(idx,:), room, i);
            intersections{idx} = intersection;
    
            if validIntersection && receiverCanSeePlane(i)
                % Check if path is blocked
                obstruction = CheckForObstruction(intersection, receiver, room, i);
                obstruction = CheckForObstruction(source, intersection, room, i, obstruction);
                if obstruction
                    validPath(idx) = false;
                else
                    count = count + 1;
                    validPath(idx) = true;
                    pathIdx(count) = idx;
                    specularPathLength(count) = pathLength;
                    refOrderPaths(count) = 1;
                end
            end
        end
        % sp-ed
        for j = 1:room.numEdges
            if edgeCanSeePlane(j,i)
                vSourceCanSeeEdge{i}(j) = ~CheckForObstruction(vSources(idx,:), room.edgeMidPoint(j,:), room, i);
                
                if receiverCanSeeEdge(j) && vSourceCanSeeEdge{i}(j)
                    % Find intersection point and check it is valid
                    [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.edgeMidPoint(j,:), vSources(idx,:), room, i);
                    %intersections{idx} = intersection;
                    
                    if validIntersection
                        obstruction = CheckForObstruction(intersection, room.edgeMidPoint(j,:), room, i);
                        obstruction = CheckForObstruction(source, intersection, room, i, obstruction);
                
                        if ~obstruction
                            diffIdx = diffIdx + 1;
                            validDiffPath(diffIdx) = true;
                            diffPath(diffIdx,:) = [j, zeros(1, diffOrder - 1)];
                            edgePath{diffIdx} = [source; intersection; room.edgeMidPoint(j,:); receiver];
                        end
                    end
                end
            end
        end
        % ed-sp
        % Check if receiver behind the plane
        [k, validReceiver] = PointPlanePosition(receiver, room.planeNormals(i,:), room.d(i));
        if validReceiver
            vR = receiver - 2 * room.planeNormals(i,:) * k;
            for j = 1:room.numEdges
                vReceiverCanSeeEdge{i}(j) = ~CheckForObstruction(vR, room.edgeMidPoint(j,:), room, i);
                if sourceCanSeeEdge(j) && vReceiverCanSeeEdge{i}(j)
                    % Find intersection point and check it is valid
                    [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.edgeMidPoint(j,:), vR, room, i);
                    %intersections{idx} = intersection;
                    
                    if validIntersection
                        obstruction = CheckForObstruction(room.edgeMidPoint(j,:), intersection, room, i);
                        obstruction = CheckForObstruction(intersection, receiver, room, i, obstruction);
                
                        if ~obstruction
                            diffIdx = diffIdx + 1;
                            validDiffPath(diffIdx) = true;
                            diffPath(diffIdx,:) = [j, zeros(1, diffOrder - 1)];
                            edgePath{diffIdx} = [source; room.edgeMidPoint(j,:); intersection; receiver];
                        end
                    end
                end
            end
        end
    end
end
toc

numVSources(1) = idx;
PlotGeometry(room.corners, room.planeCorners, source, receiver, lineOfSight, vSources, intersections, validPath, edgePath)

%% Higher order reflection

for j = 2:refOrder
    tic
    disp(['Reflection order: ', num2str(j)])
    for i = 1:room.numPlanes
        disp(['Plane: ', num2str(i)])
        for n = 1:numVSources(j - 1)
            vSourceIdx = sum(numVSources(1:j - 2)) + n;
            lastRefPlane = refPlanePath(vSourceIdx,j - 1);
            if lastRefPlane ~= i && planeCanSeeplane(i, lastRefPlane)

                vSource = vSources(vSourceIdx,:);
                % Check if source behind the plane
                [k, validSource] = PointPlanePosition(vSource, room.planeNormals(i,:), room.d(i));
                        
                if validSource
                    % Create virtual source
                    vSourceStore = vSource - 2 * room.planeNormals(i,:) * k;
                    pathLength = norm(vSourceStore - receiver);
                    if pathLength > maxPathLength
                        validSource = false;
                    end
                    if validSource
                        idx = idx + 1;
                        vSources(idx,:) = vSourceStore;

                        refPlanePath(idx,:) = [refPlanePath(vSourceIdx,1:j - 1), i, zeros(1, refOrder - j)];
                        vSourceIdxPath(idx,:) = [vSourceIdxPath(vSourceIdx,1:j - 1), idx, zeros(1, refOrder - j)];
    
                        % Find intersection point and check it is valid
                        [intersection, validPath(idx,:)] = CheckValidLinePlaneIntersection(receiver, vSources(idx,:), room, i);
                        intersections{idx}(j,:) = intersection;
                        validIntersection = validPath(idx,:);
                        m = 1;
                        while validIntersection && m < j
                            refPlane = refPlanePath(vSourceIdx, j - m);
                            startPos = vSources(vSourceIdxPath(vSourceIdx, j - m),:);
                            [intersection, validIntersection] = CheckValidLinePlaneIntersection(intersection, startPos, room, refPlane);
                            intersections{idx}(j - m,:) = intersection;
                            m = m + 1;
                        end
                        if validIntersection && receiverCanSeePlane(i)
                            % Check if path is blocked
                            obstruction = CheckForObstruction(intersections{idx}(j,:), receiver, room, i);
                            obstruction = CheckForObstruction(source, intersections{idx}(1,:), room, i, obstruction);
                            p = 1;
                            while ~obstruction && p < j
                                obstruction = CheckForObstruction(intersections{idx}(p,:), intersections{idx}(p + 1,:), room, i, obstruction);
                                p = p + 1;
                            end
                            if obstruction
                                validPath(idx) = false;
                            else
                                count = count + 1;
                                validPath(idx) = true;
                                pathIdx(count) = idx;
                                specularPathLength(count) = norm(vSources(idx,:) - receiver);
                                refOrderPaths(count) = j;
                            end
                        else
                            validPath(idx) = false;
                        end
                    end
                    % sp-ed
                    for k = 1:room.numEdges
                        if edgeCanSeePlane(k,i)
                            vSourceCanSeeEdge{i}(k) = ~CheckForObstruction(vSources(idx,:), room.edgeMidPoint(k,:), room, i);
                            
                            if receiverCanSeeEdge(k) && vSourceCanSeeEdge{i}(k)
                                % Find intersection point and check it is valid
                                [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.edgeMidPoint(k,:), vSources(idx,:), room, i);
                                intersectionStore(j,:) = intersection;
                                %intersections{idx} = intersection;
                                m = 1;
                                while validIntersection && m < j
                                    refPlane = refPlanePath(vSourceIdx, j - m);
                                    startPos = vSources(vSourceIdxPath(vSourceIdx, j - m),:);
                                    [intersection, validIntersection] = CheckValidLinePlaneIntersection(intersection, startPos, room, refPlane);
                                    intersectionStore(j - m,:) = intersection;
                                    m = m + 1;
                                end
                                if validIntersection && receiverCanSeePlane(i)
                                    % Check if path is blocked
                                    obstruction = CheckForObstruction(intersectionStore(j,:), receiver, room, i);
                                    obstruction = CheckForObstruction(source, intersectionStore(1,:), room, i, obstruction);
                                    p = 1;
                                    while ~obstruction && p < j
                                        obstruction = CheckForObstruction(intersectionStore(p,:), intersectionStore(p + 1,:), room, i, obstruction);
                                        p = p + 1;
                                    end
                                    if ~obstruction
                                        diffIdx = diffIdx + 1;
                                        validDiffPath(diffIdx) = true;
                                        diffPath(diffIdx,:) = [k, zeros(1, diffOrder - 1)];
                                        edgePath{diffIdx} = [source; intersectionStore; room.edgeMidPoint(k,:); receiver];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    toc

    PlotGeometry(room.corners, room.planeCorners, source, receiver, lineOfSight, vSources, intersections, validPath, edgePath)
    numVSources(j) = idx - sum(numVSources(1:j - 1));
end

% vSources = vSources(1:idx,:);
% validPath = validPath(1:idx);

PlotGeometry(room.corners, room.planeCorners, source, receiver, lineOfSight, vSources, intersections, validPath, edgePath)

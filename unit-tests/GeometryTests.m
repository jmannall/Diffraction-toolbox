
close all
clear all
set(0, 'DefaultLineLineWidth', 1.0);

%% Shoebox room

x = 4.97;
y = 4.12;
z = 3;

roomFunc = @()CreateShoeboxRoomGeometry(x, y, z);
room.receiver = [2.77, 1.3, 1.51];
room.source = [1.4, 3.02, 1.36];

%% Ring shaped room

% x = [8 2 5];
% y = [12 3 6];
% z = 2.6;
% 
% roomFunc = @()CreateRingShapedRoomGeometry(x, y, z);

%% L shaped room

% x = [4 7];
% y = [3 5];
% z = 3;
% 
% roomFunc = @()CreateLShapedRoomGeometry(x, y, z);
% room.receiver = [5, 4.5, 1];
% room.source = [1.5, 2. 1];

%% Modelling

refOrder = 3;
diffOrder = 3;
combinedOrder = 3;
maxPathLength = 100;

op = struct('geometry', true, 'sr', true, 'direct', false, 'specular', false, 'hodSpecular', false, 'diff', false, 'hodDiff', false);

[room, plot] = CreateRoomGeometry(roomFunc, room);

%% Direct sound

[direct, plot] = DirectComponent(room, plot);
op.direct = true;

PlotIS(plot, op)

%% First order diffraction

[diff, plot] = DiffComponenet(room, plot);
op.diff = true;

PlotIS(plot, op)

% PlotGeometry(room.corners, room.planeCorners, source, receiver, lineOfSight, vSources, intersections, validPath, edgePath)

%% Higher order diffraction

[hodDiff, plot] = HodDiffComponent(room, plot, diffOrder);
op.hodDiff = true;

PlotIS(plot, op)

%% First order reflections

[specular, plot, sData] = SpecularComponent(room, plot);
op.specular = true;

PlotIS(plot, op)

%% Higher order reflections

[hodSpecular, plot, sData] = HodSpecularComponent(room, plot, sData, refOrder);
op.hodSpecular = true;

PlotIS(plot, op)





%vSources = zeros(inf, 3);
%validPath = false(inf, 1);
%[pathIdx, specularPathLength, refOrderPaths] = deal(zeros(inf, 1));
%[refPlanePath, vSourceIdxPath] = deal(zeros(inf, refOrder));

diffIdx = 0;
idx = 0;
count = 0;
tic
disp('Reflection order: 1')
for i = 1:room.numPlanes
    disp(['Plane: ', num2str(i)])
    % Check if source behind the plane
    [k, validSource] = PointPlanePosition(room.source, room.planeNormals(i,:), room.d(i));

    if validSource
        % Create virtual source
        vS = room.source - 2 * room.planeNormals(i,:) * k;
        pathLength = PathLength(vS, room.receiver);
        if pathLength > maxPathLength
            validSource = false;
        end
        if validSource
            idx = idx + 1;  % Source will be stored
            vSources(idx,:) = vS;
            refPlanePath(idx,:) = [i, zeros(1, refOrder - 1)];
            vSourceIdxPath(idx,:) = [idx, zeros(1, refOrder - 1)];
    
            % Find intersection point and check it is valid
            [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.receiver, vSources(idx,:), room, i);
            intersections{idx} = intersection;
    
            if validIntersection && room.receiverCanSeePlane(i)
                % Check if path is blocked
                obstruction = CheckForObstruction(intersection, room.receiver, room, i);
                obstruction = CheckForObstruction(room.source, intersection, room, i, obstruction);
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
            if room.edgeCanSeePlane(j,i)
                vSourceCanSeeEdge{i}(j) = ~CheckForObstruction(vSources(idx,:), room.edgeMidPoints(j,:), room, i);
                
                if room.receiverCanSeeEdge(j) && vSourceCanSeeEdge{i}(j)
                    % Find intersection point and check it is valid
                    [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.edgeMidPoints(j,:), vSources(idx,:), room, i);
                    %intersections{idx} = intersection;
                    
                    if validIntersection
                        obstruction = CheckForObstruction(intersection, room.edgeMidPoints(j,:), room, i);
                        obstruction = CheckForObstruction(room.source, intersection, room, i, obstruction);
                
                        if ~obstruction
                            diffIdx = diffIdx + 1;
                            validDiffPath(diffIdx) = true;
                            diffPath(diffIdx,:) = [j, zeros(1, diffOrder - 1)];
                            edgePath{diffIdx} = [room.source; intersection; room.edgeMidPoints(j,:); room.receiver];
                        end
                    end
                end
            end
        end
        % ed-sp
        % Check if receiver behind the plane
        [k, validReceiver] = PointPlanePosition(room.receiver, room.planeNormals(i,:), room.d(i));
        if validReceiver
            vR = room.receiver - 2 * room.planeNormals(i,:) * k;
            for j = 1:room.numEdges
                vReceiverCanSeeEdge{i}(j) = ~CheckForObstruction(vR, room.edgeMidPoints(j,:), room, i);
                if room.sourceCanSeeEdge(j) && vReceiverCanSeeEdge{i}(j)
                    % Find intersection point and check it is valid
                    [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.edgeMidPoints(j,:), vR, room, i);
                    %intersections{idx} = intersection;
                    
                    if validIntersection
                        obstruction = CheckForObstruction(room.edgeMidPoints(j,:), intersection, room, i);
                        obstruction = CheckForObstruction(intersection, room.receiver, room, i, obstruction);
                
                        if ~obstruction
                            diffIdx = diffIdx + 1;
                            validDiffPath(diffIdx) = true;
                            diffPath(diffIdx,:) = [j, zeros(1, diffOrder - 1)];
                            edgePath{diffIdx} = [room.source; room.edgeMidPoints(j,:); intersection; room.receiver];
                        end
                    end
                end
            end
        end
    end
end
toc

numVSources(1) = idx;
PlotGeometry(room.corners, room.planeCorners, room.source, room.receiver, lineOfSight, vSources, intersections, validPath, edgePath)

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
                            vSourceCanSeeEdge{i}(k) = ~CheckForObstruction(vSources(idx,:), room.edgeMidPoints(k,:), room, i);
                            
                            if receiverCanSeeEdge(k) && vSourceCanSeeEdge{i}(k)
                                % Find intersection point and check it is valid
                                [intersection, validIntersection] = CheckValidLinePlaneIntersection(room.edgeMidPoints(k,:), vSources(idx,:), room, i);
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
                                        edgePath{diffIdx} = [source; intersectionStore; room.edgeMidPoints(k,:); receiver];
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

close all
clear all

x = [2 4];
y = [3 5];
z = 2.5;

[corners, planes, source, receiver] = CreateLShapedRoomGeometry(x, y, z);

source = [1.1 0.9 1];
receiver = [2.4 3.8 1];

refOrder = 6;

% corners = [0 0 0
%     x(1) 0 0
%     x(1) y(1) 0
%     0 y(1) 0
%     0 0 z
%     x(1) 0 z
%     x y(1) z
%     0 y(1) z];
% 
% edges = [1 2
%     2 3
%     3 4
%     4 1
%     5 6
%     6 7
%     7 8
%     8 5
%     1 5
%     2 6
%     3 7
%     4 8];
% 
% planes = [1 2 3 4
%     5 8 7 6
%     1 5 6 2
%     2 6 7 3
%     3 7 8 4
%     4 8 5 1];

numPlanes = size(planes, 1);
normals = zeros(numPlanes, 3);
[d, numPlaneCorners] = deal(zeros(numPlanes, 1));
validPlane = true(numPlanes, 1);

% Create geometry parameters
for i = 1:numPlanes
    validCorners = planes(i,:) > 0;
    numPlaneCorners(i) = sum(validCorners);
    planeCorners = corners(planes(i,validCorners),:);
    normals(i,:) = normr(cross(planeCorners(1,:) - planeCorners(2,:), planeCorners(1,:) - planeCorners(3,:)));
    d(i) = sum(normals(i,:) .* planeCorners(1,:));
    % Check if receiver lies behind the plane
    [kR(i), validPlane(i)] = PointPlanePosition(receiver, normals(i,:), d(i));
end
edgesCanSee = true(numPlanes);
for i = 1:numPlanes
    for j = 1:numPlanes
        dotProduct = dot(normals(i,:), normals(j,:));
        angle1 = acosd(dot(normals(i,:), normals(j,:)));
        angle2 = asind(dot(normals(i,:), normals(j,:)));
        if dotProduct == 1
            [edgesCanSee(i,j), edgesCanSee(j,i)] = deal(false);
        else
            for n = 1:numPlaneCorners(j)
                validCorners = false(1, numPlaneCorners(j));
                [~, validCorners(n)] = PointPlanePosition(corners(planes(j,n),:), normals(i,:), d(i));
            end
            if sum(validCorners) == 0
                [edgesCanSee(i,j), edgesCanSee(j,i)] = deal(false);
            end
        end
    end
end

% First order reflections
idx = 0;
count = 0;
for i = 1:numPlanes
    % Check if source behind the plane
    [k, validSource] = PointPlanePosition(source, normals(i,:), d(i));

    if validSource
        % Create virtual source
        idx = idx + 1;
        vSources(idx,:) = source - 2 * normals(i,:) * k;
        refPlanePath(idx,:) = [i, zeros(1, refOrder - 1)];
        vSourceIdxPath(idx,:) = [idx, zeros(1, refOrder - 1)];

        % Find intersection point and check it is valid
        [intersection, validPath(idx,:)] = CheckValidLinePlaneIntersection(receiver, vSources(idx,:), normals(i,:), d(i), planes, corners, numPlaneCorners, i);
        intersections{idx} = intersection;

        if validPath(idx)
            if validPlane(i)
                % Check if path is blocked
                obstruction = false;
                obstruction = CheckForObstruction(intersection, receiver, planes, normals, d, corners, numPlaneCorners, numPlanes, i, obstruction);
                obstruction = CheckForObstruction(source, intersection, planes, normals, d, corners, numPlaneCorners, numPlanes, i, obstruction);
                if obstruction
                    validPath(idx) = false;
                else
                    count = count + 1;
                    validPath(idx) = true;
                    specularPathLength(count) = norm(vSources(idx,:) - receiver);
                    refOrderPaths(count) = 1;
                end
            end
        end
    end
end
PlotGeometry(corners, planes, source, receiver, vSources, intersections, validPath)

% Correct up to first order reflection not including accounting for
% obstructions
% Higher order reflections
numVSources(1) = idx;
for j = 2:refOrder
    for i = 1:numPlanes
        for n = 1:numVSources(j - 1)
            vSourceIdx = sum(numVSources(1:j - 2)) + n;
            lastRefPlane = refPlanePath(vSourceIdx,j - 1);
            if lastRefPlane ~= i && edgesCanSee(i, lastRefPlane)

                vSource = vSources(vSourceIdx,:);
                % Check if source behind the plane
                [k, validSource] = PointPlanePosition(vSource, normals(i,:), d(i));
                        
                if validSource
                    % Create virtual source
                    idx = idx + 1;
                    vSources(idx,:) = vSource - 2 * normals(i,:) * k;
                    refPlanePath(idx,:) = [refPlanePath(vSourceIdx,1:j - 1), i, zeros(1, refOrder - j)];
                    vSourceIdxPath(idx,:) = [vSourceIdxPath(vSourceIdx,1:j - 1), idx, zeros(1, refOrder - j)];

                    % Find intersection point and check it is valid
                    [intersection, validPath(idx,:)] = CheckValidLinePlaneIntersection(receiver, vSources(idx,:), normals(i,:), d(i), planes, corners, numPlaneCorners, i);
                    intersections{idx}(j,:) = intersection;
                    for m = 1:j - 1
                        refPlane = refPlanePath(vSourceIdx, j - m);
                        startPos = vSources(vSourceIdxPath(vSourceIdx, j - m),:);
                        [intersection, validPathStore] = CheckValidLinePlaneIntersection(intersection, startPos, normals(refPlane,:), d(refPlane), planes, corners, numPlaneCorners, refPlane);
                        intersections{idx}(j - m,:) = intersection;
                        if ~validPathStore
                            validPath(idx) = validPathStore;
                        end
                    end

                    if validPath(idx)
                        if validPlane(i)
                            % Check if path is blocked
                            obstruction = false;
                            obstruction = CheckForObstruction(intersections{idx}(j,:), receiver, planes, normals, d, corners, numPlaneCorners, numPlanes, i, obstruction);
                            obstruction = CheckForObstruction(source, intersections{idx}(1,:), planes, normals, d, corners, numPlaneCorners, numPlanes, i, obstruction);
                            p = 1;
                            while ~obstruction && p < j
                                obstruction = CheckForObstruction(intersections{idx}(p,:), intersections{idx}(p + 1,:), planes, normals, d, corners, numPlaneCorners, numPlanes, i, obstruction);
                                p = p + 1;
                            end
                            if obstruction
                                validPath(idx) = false;
                            else
                                count = count + 1;
                                validPath(idx) = true;
                                specularPathLength(count) = norm(vSources(idx,:) - receiver);
                                refOrderPaths(count) = j;
                            end
                        else
                            validPath(idx) = false;
                        end
                    end
                end
            end
        end
    end
    PlotGeometry(corners, planes, source, receiver, vSources, intersections, validPath)
    numVSources(j) = idx - sum(numVSources(1:j - 1));
end

numValidPaths = sum(validPath);
refPlanePathsValid = refPlanePath(validPath,:);
% Equation of a plane is ax + by + cz = d where a, b, c are normal(x, y, z)
% and d is sum(normal(x, y, z) .* pointOnPlane(x, y, z)).

%% Audio

c = 344;
fs = 96e3;
nfft = 8192;
windowLength = 100;
validPath = true;

audio = audioread('audio/Imagination.wav');
tfcomplexAll = zeros(nfft / 2, 1);
for i = 1:length(specularPathLength)
    [output, ir] = DelayLine(audio, specularPathLength(i), windowLength, validPath, c, fs);
    [~, tfcomplexStore] = IrToTf(ir, nfft);
    tfcomplexAll = tfcomplexAll + tfcomplexStore;
end

tfcomplexThree = zeros(nfft / 2, 1);
specularPathLengthThree = specularPathLength(refOrderPaths < 4);
for i = 1:length(specularPathLengthThree)
    [output, ir] = DelayLine(audio, specularPathLengthThree(i), windowLength, validPath, c, fs);
    [~, tfcomplexStore] = IrToTf(ir, nfft);
    tfcomplexThree = tfcomplexThree + tfcomplexStore;
end

fvec = fs/nfft*[0:nfft/2-1];
tfmagAll = mag2db(abs(tfcomplexAll));
tfmagThree = mag2db(abs(tfcomplexThree));

figure
semilogx(fvec, tfmagAll)
hold on
semilogx(fvec, tfmagThree)
legend('All', 'Three')
xlim([20 20e3])
ylim([-70 20])
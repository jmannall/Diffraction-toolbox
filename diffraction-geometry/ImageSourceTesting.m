close all
clear all

file = 'audio\ESS.wav';
fs = 44.1e3;
duration = 1.5;
scilence = 2;
ess = sweeptone(duration,scilence,fs);

audiowrite(file, ess, fs)

x = [2 4];
y = [3 5];
z = 2.5;

[corners, planes, source, receiver] = CreateLShapedRoomGeometry(x, y, z);

x = [20 1.5 18.5];
y = [24 2 19];
z = 3;

receiver = [19.25, 1, 1.5];
[corners, planes, source, receiver] = CreateRingShapedRoomGeometry(x, y, z, receiver);

%source = [1.1 0.9 1];
%receiver = [3 3.8 1];

refOrder = 6;
maxPathLength = 100;

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
            validCorners = false(1, numPlaneCorners(j));
            for n = 1:numPlaneCorners(j)
                [~, validCorners(n)] = PointPlanePosition(corners(planes(j,n),:), normals(i,:), d(i));
            end
            if sum(validCorners) == 0
                [edgesCanSee(i,j), edgesCanSee(j,i)] = deal(false);
            end
        end
    end
end


% Check for direct sound

i = 0;
obstruction = false;
obstruction = CheckForObstruction(source, receiver, planes, normals, d, corners, numPlaneCorners, numPlanes, i, obstruction);

lineOfSight = ~obstruction;
if lineOfSight
    directPathLength = PathLength(receiver, source);
end

% First order reflections
idx = 0;
count = 0;
disp('Reflection order: 1')
for i = 1:numPlanes
    disp(['Plane: ', num2str(i)])
    % Check if source behind the plane
    [k, validSource] = PointPlanePosition(source, normals(i,:), d(i));

    if validSource
        % Create virtual source
        vSourceStore = source - 2 * normals(i,:) * k;
        pathLength = norm(vSourceStore - receiver);
        if pathLength > maxPathLength
            validSource = false;
        end
        if validSource
            idx = idx + 1;
            vSources(idx,:) = vSourceStore;
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
                        pathIdx(count) = idx;
                        specularPathLength(count) = norm(vSources(idx,:) - receiver);
                        refOrderPaths(count) = 1;
                    end
                end
            end
        end
    end
end
PlotGeometry(corners, planes, source, receiver, lineOfSight, vSources, intersections, validPath)

% Correct up to first order reflection not including accounting for
% obstructions
% Higher order reflections
numVSources(1) = idx;
for j = 2:refOrder
    disp(['Reflection order: ', num2str(j)])
    for i = 1:numPlanes
        disp(['Plane: ', num2str(i)])
        for n = 1:numVSources(j - 1)
            vSourceIdx = sum(numVSources(1:j - 2)) + n;
            lastRefPlane = refPlanePath(vSourceIdx,j - 1);
            if lastRefPlane ~= i && edgesCanSee(i, lastRefPlane)

                vSource = vSources(vSourceIdx,:);
                % Check if source behind the plane
                [k, validSource] = PointPlanePosition(vSource, normals(i,:), d(i));
                        
                if validSource
                    % Create virtual source
                    vSourceStore = vSource - 2 * normals(i,:) * k;
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
                                    pathIdx(count) = idx;
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
    end
    PlotGeometry(corners, planes, source, receiver, lineOfSight, vSources, intersections, validPath)
    numVSources(j) = idx - sum(numVSources(1:j - 1));
end

numValidPaths = sum(validPath);
refPlanePathsValid = refPlanePath(validPath,:);
% Equation of a plane is ax + by + cz = d where a, b, c are normal(x, y, z)
% and d is sum(normal(x, y, z) .* pointOnPlane(x, y, z)).

%% Audio

close all

PlotGeometry(corners, planes, source, receiver, lineOfSight, vSources, intersections, validPath)

audioFile = 'piano';
[audio, fs] = audioread(['sourceAudio/' audioFile '.wav']);

scale = 20;
audio = audio * scale;
c = 344;
% fs = 96e3;
nfft = 8192;
windowLength = 100;
controlparameters = struct('fs', fs, 'c', c, 'nfft', nfft, 'difforder', 2, 'saveFiles', 1, 'noDirect', false);

heading = [0 -1 0];
planeCorners = [1 4 3 2
    5 6 7 8
    1 2 6 5
    2 3 7 6
    3 4 8 7
    4 1 5 8];
planeRigid = [1 1 1 1 0 1];

[ir, tfmag, tvec, fvec, tfcomplex] = SingleBTM(source, receiver, corners([5:8, 13:16],:), planeCorners, planeRigid, controlparameters, true);
[azimuthD(1,1), elevationD(1,1)] = CalculateAzimuthElevation(heading, receiver, [x(2) y(2) z / 2]);

irD{1} = ir.diff2;

planeRigid = [1 1 0 1 1 1];

[ir, tfmag, tvec, fvec, tfcomplex] = SingleBTM(source, receiver, corners([5:8, 13:16],:), planeCorners, planeRigid, controlparameters, true);
[azimuthD(2,1), elevationD(2,1)] = CalculateAzimuthElevation(heading, receiver, [x(2) y(3) z / 2]);

irD{2} = ir.diff2;

%heading = zeros(size([source; vSources(pathIdx,:)])) .* NormaliseVector(source - receiver);
heading = ones(size([source; vSources(pathIdx,:)])) .* (heading);

[azimuth, elevation] = CalculateAzimuthElevation(heading, receiver, [source; vSources(pathIdx,:)]);

azimuth = [azimuth; azimuthD];
elevation = [elevation; elevationD];

hrtf = CreateHRTF(azimuth, elevation);

tfcomplexAll = zeros(nfft / 2, 1);
[delay, fracDelay] = CalculateDelay(specularPathLength, c, fs);
irLength = max(max(delay) + windowLength, max(length(irD{1}), length(irD{2}))) + 256;
[brirL, brirR] = deal(zeros(irLength, 1));

if lineOfSight
    start = 0;
else
    start = 1;
end

for i = start:length(specularPathLength)
    if i > 0
        pathLength = specularPathLength(i);
        alpha = 0.95;
    else
        pathLength = directPathLength;
        alpha = 1;
    end
    [output, ir] = DelayLine(audio, pathLength, windowLength, true, c, fs);
    ir = alpha * ir;
    [~, tfcomplexStore] = IrToTf(ir, nfft);
    tfcomplexAll = tfcomplexAll + tfcomplexStore;

    birL = conv(ir, hrtf.L(:,i + 1));
    birLength = length(birL);
    brirL(1:birLength) = brirL(1:birLength) + birL;

    birR = conv(ir, hrtf.R(:,i + 1));
    birLength = length(birR);
    brirR(1:birLength) = brirR(1:birLength) + birR;
end

outputL = conv(audio, brirL);
outputR = conv(audio, brirR);

audiowrite(['audio/' audioFile '_All.wav'], [outputL, outputR], fs)

brirLRev = brirL;
brirRRev = brirR;

birL = conv(irD{1}, hrtf.L(:,i + 2));
birLength = length(birL);
brirL(1:birLength) = brirL(1:birLength) + birL;

birR = conv(irD{1}, hrtf.R(:,i + 2));
birLength = length(birR);
brirR(1:birLength) = brirR(1:birLength) + birR;

birL = conv(irD{2}, hrtf.L(:,i + 3));
birLength = length(birL);
brirL(1:birLength) = brirL(1:birLength) + birL;

birR = conv(irD{2}, hrtf.R(:,i + 3));
birLength = length(birR);
brirR(1:birLength) = brirR(1:birLength) + birR;

[~, tfcomplexStore] = IrToTf(irD{1}, nfft);
tfcomplexDiff = tfcomplexAll + tfcomplexStore;
tfcomplexDiffOnly = tfcomplexStore;

[~, tfcomplexStore] = IrToTf(irD{2}, nfft);
tfcomplexDiff = tfcomplexDiff + tfcomplexStore;
tfcomplexDiffOnly = tfcomplexDiffOnly + tfcomplexStore;

outputL = conv(audio, brirL);
outputR = conv(audio, brirR);

audiowrite(['audio/' audioFile '_AllDiff.wav'], [outputL, outputR], fs)

fvec = fs/nfft*[0:nfft/2-1];
tfmagAll = mag2db(abs(tfcomplexAll));
tfmagDiff = mag2db(abs(tfcomplexDiff));
tfmagDiffOnly = mag2db(abs(tfcomplexDiffOnly));


t = (1:length(impulse_response)) / 44.1e3;
figure
plot(t, impulse_response(:,1))
grid on
title('Left rtSOFE')
xlim([0 0.15])
ylim([-0.01 0.01])

figure
plot(t, impulse_response(:,2))
grid on
title('Right rtSOFE')
xlim([0 0.15])
ylim([-0.01 0.01])

figure
semilogx(fvec, tfmagAll)
hold on
grid on
semilogx(fvec, tfmagDiff)
semilogx(fvec, tfmagDiffOnly)
legend('All', 'Diff', 'DiffOnly')
xlim([20 20e3])
ylim([-70 20])

tvec = (1:length(brirL)) / fs;
figure
% tiledlayout(2,1)

% nexttile
plot(tvec, brirLRev)
grid on
title('Left IS')
ylim([-0.01 0.01])

% nexttile
% plot(tvec, brirL)
% grid on

figure
% tiledlayout(2,1)

% nexttile
plot(tvec, brirRRev)
grid on
title('Right IS')
ylim([-0.01 0.01])

% nexttile
% plot(tvec, brirR)
% grid on

figure
tiledlayout(2,1)

nexttile
plot(tvec, brirLRev + brirRRev)
grid on

nexttile
plot(tvec, brirL + brirR)
grid on


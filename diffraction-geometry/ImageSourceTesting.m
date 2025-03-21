close all
clear all







i = 1:6;
count = [0 0 0];
idx2 = 1;
idx3 = 1;
for k = 1:length(i)
   count(1) = count(1) + 1;
   path{1}(k) = k;
   for j = 1:length(i)
        if j ~= k
            count(2) = count(2) + 1;
            path{2}(idx2,:) = [k, j];
            idx2 = idx2 + 1;
            for m = 1:length(i)
                if m ~= j
                    count(3) = count(3) + 1;
                    path{3}(idx3,:) = [k, j, m];
                    idx3 = idx3 + 1;
                end
            end
        end
   end
end

path{2} = sort(path{2}, 2);

audioFile = 'music';
[audio, fs] = audioread(['sourceAudio/' audioFile '.wav']);

x = [1 3];
y = [3 5];
z = 2;

room = CreateLShapedRoomGeometry(x, y, z);

source = [2.5, 4.4, 1];
receiver = [0.7, 1.8, 1.1];

x = [15.73 1.73 14];
y = [23.72 2.33 19.58];
z = 2.6;

% x = [4.97 2 3];
% y = [4.12 2 3];
% z = 3;
% 
% [corners, planes, source] = CreateRingShapedRoomGeometry(x, y, z);
% 
% source = [1.40, 3.02, 1.36];
% receiver = [2.77, 1.30, 1.51];

x = 3;
y = 5;
z = 2;
% room = CreateShoeboxRoomGeometry(x, y, z);
% 
% 
% source = [2.5 3.23 1];
% receiver = [1.95 2.1 1.1];

% source = [2 3 1];
% receiver = [2 2 1];
% x = [2 3];
% y = [2 4];
% z = 2; 
% room = CreateLShapedRoomGeometry(x, y, z);
room.source = source;
room.receiver = receiver;
PlotGeometry(room.corners, room.planeCorners, source, receiver)

alphaFreq = [250 500 1e3 2e3 4e3];

alpha = [0.06 0.15 0.4 0.6 0.6  % floor
    0.01 0.02 0.02 0.02 0.03    % ceiling
    0.1 0.05 0.04 0.07 0.1      % -y
    0.3 0.1 0.1 0.1 0.1         % -x
    0.7 0.6 0.7 0.7 0.5         % y
    0.2 0.2 0.1 0.07 0.04];     % x

% f1 = 20;
% f2 = max(alphaFreq)^(0.1) * (fs/2)^(0.9);
% 
% [bpB, bpA] = butter(1, [f1 ,f2]/fs*2);
% 
% refBTemp = ism_setup.b_refl;
% refATemp = ism_setup.a_refl;
% 
% idx = [1 6 2 3 5 4];
% for i = 1:6
%     refB{1,i} = refBTemp{1,idx(i)};
%     refA{1,i} = refATemp{1,idx(i)};
% 
%     refB{2,i} = refBTemp{2,idx(i)};
%     refA{2,i} = refATemp{2,idx(i)};
% end

%source = [1.1 0.9 1];
%receiver = [3 3.8 1];
tic
refOrder = 3;
maxPathLength = 1e3;

room.numPlanes = size(room.planeCorners, 1);
room.planeNormals = zeros(room.numPlanes, 3);
[room.d, room.numPlaneCorners] = deal(zeros(room.numPlanes, 1));
validPlane = true(room.numPlanes, 1);

% Create geometry parameters
for i = 1:room.numPlanes
    validCorners = room.planeCorners(i,:) > 0;
    room.numPlaneCorners(i) = sum(validCorners);
    planeCorners = room.corners(room.planeCorners(i,validCorners),:);
    room.planeNormals(i,:) = normr(cross(planeCorners(1,:) - planeCorners(2,:), planeCorners(1,:) - planeCorners(3,:)));
    room.d(i) = sum(room.planeNormals(i,:) .* planeCorners(1,:));
    % Check if receiver lies behind the plane
    [kR(i), validPlane(i)] = PointPlanePosition(receiver, room.planeNormals(i,:), room.d(i));
end
edgesCanSee = true(room.numPlanes);
for i = 1:room.numPlanes
    for j = 1:room.numPlanes
        dotProduct = dot(room.planeNormals(i,:), room.planeNormals(j,:));
        angle1 = acosd(dot(room.planeNormals(i,:), room.planeNormals(j,:)));
        angle2 = asind(dot(room.planeNormals(i,:), room.planeNormals(j,:)));
        if dotProduct == 1
            [edgesCanSee(i,j), edgesCanSee(j,i)] = deal(false);
        else
            validCorners = false(1, room.numPlaneCorners(j));
            for n = 1:room.numPlaneCorners(j)
                [~, validCorners(n)] = PointPlanePosition(room.corners(room.planeCorners(j,n),:), room.planeNormals(i,:), room.d(i));
            end
            if sum(validCorners) == 0
                [edgesCanSee(i,j), edgesCanSee(j,i)] = deal(false);
            end
        end
    end
end
%%

% Check for direct sound
i = 0;
obstruction = false;
obstruction = CheckForObstruction(source, receiver, room, i, obstruction);

lineOfSight = ~obstruction;
if lineOfSight
    directPathLength = PathLength(receiver, source);
end

% First order reflections
idx = 0;
count = 0;
disp('Reflection order: 1')
for i = 1:room.numPlanes
    disp(['Plane: ', num2str(i)])
    % Check if source behind the plane
    [k, validSource] = PointPlanePosition(room.source, room.planeNormals(i,:), room.d(i));

    if validSource
        % Create virtual source
        vSourceStore = room.source - 2 * room.planeNormals(i,:) * k;
        pathLength = norm(vSourceStore - room.receiver);
%         if pathLength > maxPathLength
%             validSource = false;
%         end
        if validSource
            idx = idx + 1;
            vSources(idx,:) = vSourceStore;
            refPlanePath(idx,:) = [i, zeros(1, refOrder - 1)];
            vSourceIdxPath(idx,:) = [idx, zeros(1, refOrder - 1)];
    
            % Find intersection point and check it is valid
            [intersection, validPath(idx,:)] = CheckValidLinePlaneIntersection(receiver, vSources(idx,:), room, i);
            intersections{idx} = intersection;
    
            if validPath(idx)
                if validPlane(i)
                    % Check if path is blocked
                    obstruction = false;
                    obstruction = CheckForObstruction(intersection, receiver, room, i, obstruction);
                    obstruction = CheckForObstruction(source, intersection, room, i, obstruction);
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
toc
PlotGeometry(room.corners, room.planeCorners, room.source, room.receiver, lineOfSight, vSources, intersections, validPath)

% Correct up to first order reflection not including accounting for
% obstructions
% Higher order reflections
numVSources(1) = idx;
numVisVSources(1) = count;
for j = 2:refOrder
    disp(['Reflection order: ', num2str(j)])
    for i = 1:room.numPlanes
        disp(['Plane: ', num2str(i)])
        for n = 1:numVSources(j - 1)
            vSourceIdx = sum(numVSources(1:j - 2)) + n;
            lastRefPlane = refPlanePath(vSourceIdx,j - 1);
            %if lastRefPlane ~= i && edgesCanSee(i, lastRefPlane)
            if lastRefPlane ~= i

                vSource = vSources(vSourceIdx,:);
                % Check if source behind the plane
                [k, validSource] = PointPlanePosition(vSource, room.planeNormals(i,:), room.d(i));
                        
                if validSource
                    % Create virtual source
                    vSourceStore = vSource - 2 * room.planeNormals(i,:) * k;
                    pathLength = norm(vSourceStore - room.receiver);
%                     if pathLength > maxPathLength
%                         validSource = false;
%                     end
                    if validSource
                        idx = idx + 1;
                        vSources(idx,:) = vSourceStore;

                        refPlanePath(idx,:) = [refPlanePath(vSourceIdx,1:j - 1), i, zeros(1, refOrder - j)];
                        vSourceIdxPath(idx,:) = [vSourceIdxPath(vSourceIdx,1:j - 1), idx, zeros(1, refOrder - j)];
    
                        % Find intersection point and check it is valid
                        [intersection, validPath(idx,:)] = CheckValidLinePlaneIntersection(room.receiver, vSources(idx,:), room, i);
                        intersections{idx}(j,:) = intersection;
                        for m = 1:j - 1
                            refPlane = refPlanePath(vSourceIdx, j - m);
                            startPos = vSources(vSourceIdxPath(vSourceIdx, j - m),:);
                            [intersection, validPathStore] = CheckValidLinePlaneIntersection(intersection, startPos, room, refPlane);
                            intersections{idx}(j - m,:) = intersection;
                            if ~validPathStore
                                validPath(idx) = validPathStore;
                            end
                        end
    
                        if validPath(idx)
                            if validPlane(i)
                                % Check if path is blocked
                                obstruction = false;
                                obstruction = CheckForObstruction(intersections{idx}(j,:), room.receiver, room, i, obstruction);
                                refPlane = refPlanePath(vSourceIdx, 1);
                                obstruction = CheckForObstruction(room.source, intersections{idx}(1,:), room, refPlane, obstruction);
                                p = 1;
                                while ~obstruction && p < j
                                    refPlanes = refPlanePath(idx, [j - p, j - p + 1]);
                                    obstruction = CheckForObstruction(intersections{idx}(p,:), intersections{idx}(p + 1,:), room, refPlanes, obstruction);
                                    p = p + 1;
                                end
                                if obstruction
                                    validPath(idx) = false;
                                else
                                    count = count + 1;
                                    validPath(idx) = true;
                                    pathIdx(count) = idx;
                                    specularPathLength(count) = norm(vSources(idx,:) - room.receiver);
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
    PlotGeometry(room.corners, room.planeCorners, room.source, room.receiver, lineOfSight, vSources, intersections, validPath)
    numVSources(j) = idx - sum(numVSources(1:j - 1));
    numVisVSources(j) = count - sum(numVisVSources(1:j - 1));
end

numValidPaths = sum(validPath);
refPlanePathsValid = refPlanePath(validPath,:);
% Equation of a plane is ax + by + cz = d where a, b, c are normal(x, y, z)
% and d is sum(normal(x, y, z) .* pointOnPlane(x, y, z)).

%% Audio

close all

PlotGeometry(corners, planes, source, receiver, lineOfSight, vSources, intersections, validPath)
scale = 1;
audio = audio * scale;
c = 344;
% fs = 96e3;
nfft = 8192;
windowLength = 100;
controlparameters = struct('fs', fs, 'c', c, 'nfft', nfft, 'difforder', 2, 'saveFiles', 1, 'noDirect', false);

heading = [0 1 0];
% planeCorners = [1 4 3 2
%     5 6 7 8
%     1 2 6 5
%     2 3 7 6
%     3 4 8 7
%     4 1 5 8];
% planeRigid = [0 0 1 1 1 1];
% 
% [ir, tfmag, tvec, fvec, tfcomplex] = SingleBTM(source, receiver, corners([5:8, 13:16],:), planeCorners, planeRigid, controlparameters, true);
% [azimuthD(1,1), elevationD(1,1)] = CalculateAzimuthElevation(heading, receiver, [x(2) y(2) z / 2]);
% [azimuthD(2,1), elevationD(2,1)] = CalculateAzimuthElevation(heading, receiver, [x(2) y(3) z / 2]);
% 
% irD{1} = ir.diff1;
% irD{2} = ir.diff2;

%planeRigid = [0 0 0 1 1 1];

%[ir, tfmag, tvec, fvec, tfcomplex] = SingleBTM(source, receiver, corners([5:8, 13:16],:), planeCorners, planeRigid, controlparameters, true);
%[azimuthD(2,1), elevationD(2,1)] = CalculateAzimuthElevation(heading, receiver, [x(2) y(3) z / 2]);

%irD{2} = ir.diff2;

%heading = zeros(size([source; vSources(pathIdx,:)])) .* NormaliseVector(source - receiver);
heading = ones(size([source; vSources(pathIdx,:)])) .* (heading);

[azimuth, elevation] = CalculateAzimuthElevation(heading, receiver, [source; vSources(pathIdx,:)]);

% azimuth = [azimuth; azimuthD];
% elevation = [elevation; elevationD];

hrtf = CreateHRTF(azimuth, elevation);

tfcomplexAll = zeros(nfft / 2, 1);
[delay, fracDelay] = CalculateDelay(specularPathLength, c, fs);
%irLength = max(max(delay) + windowLength, max(length(irD{1}), length(irD{2}))) + 256;
irLength = max(delay) + windowLength + 256 + 1e4;

[brirL, brirR, irStore] = deal(zeros(irLength, 1));

if lineOfSight
    start = 0;
else
    start = 1;
end

for i = start:length(specularPathLength)
    if i > 0
        pathLength = specularPathLength(i);
    else
        pathLength = directPathLength;
    end
    delay = pathLength * fs / c;
    pathLength = round(delay) * c / fs;
    [output, ir] = DelayLine(audio, pathLength, windowLength, true, c, fs);
    if i > 0
        ir = [ir; zeros(1e4, 1)];
        numRef = find(refPlanePathsValid(i,:), 1, "last");
        for j = 1:numRef
            ir = filter(refB{1,refPlanePathsValid(i,j)}, refA{1,refPlanePathsValid(i,j)}, ir);
            ir = filter(refB{2,refPlanePathsValid(i,j)}, refA{2,refPlanePathsValid(i,j)}, ir);
        end
    end
    [~, tfcomplexStore] = IrToTf(ir, nfft);
    tfcomplexAll = tfcomplexAll + tfcomplexStore;

    len = length(ir);
    irStore(1:len) = irStore(1:len) + ir;

    irAll{i+1} = ir;

    birL = conv(ir, hrtf.L(:,i + 1));
    birLength = length(birL);
    brirL(1:birLength) = brirL(1:birLength) + birL;

    birR = conv(ir, hrtf.R(:,i + 1));
    birLength = length(birR);
    brirR(1:birLength) = brirR(1:birLength) + birR;
end

irStoreBP = filter(bpB, bpA, irStore);

%irStoreBP = DelayLineIIRFilter(irStore, 0, length(irStore), true, bpB, bpA, c, fs, true);

outputL = conv(audio, brirL);
outputR = conv(audio, brirR);

audiowrite(['audio/' audioFile '_Laboratory.wav'], [outputL, outputR], fs)

tvec = (1:length(brirL)) / fs;
figure
plot(tvec, brirL)
hold on
grid on
plot(tvec, brirR)
legend('Left', 'Right')
title('Laboratoy ISM')

tvec = (1:length(irStore)) / fs;
figure
plot(tvec, irStore)
grid on
legend('ir')
title('Laboratoy ISM (Not spatialised)')
xlim([0 0.2])
ylim([-0.1 0.5])

figure
% tvec = (1:length(irD{1})) / fs;
% plot(tvec, irD{1})
hold on
% tvec = (1:length(irD{2})) / fs;
% plot(tvec, irD{2})
for i = 1:length(irAll)
    tvec = (1:length(irAll{i})) / fs;
    plot(tvec, irAll{i})
    grid on
    %legend('ir')
    title('Laboratoy ISM (Not spatialised)')
    xlim([0 0.2])
    ylim([-0.05 0.45])
end

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

% t = (1:length(impulse_response)) / 44.1e3;
% figure
% plot(t, impulse_response(:,1))
% grid on
% title('Left rtSOFE')
% xlim([0 0.15])
% ylim([-0.01 0.01])
% 
% figure
% plot(t, impulse_response(:,2))
% grid on
% title('Right rtSOFE')
% xlim([0 0.15])
% ylim([-0.01 0.01])

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

%%
load(['RIR' filesep 'BCCorridor_1stOrder']);
irData = measurementData{"omni_ESS", "ImpulseResponse"};

figure
plot(irData.Time, irData.Amplitude)
title('Omni Measured')

irDataL = measurementData{"binaural_L_ESS", "ImpulseResponse"};
figure
plot(irDataL.Time, irDataL.Amplitude)
title('L Measured')
grid on
xlim([0.1 0.2])
ylim([-0.04 0.04])

irDataR = measurementData{"binaural_R_ESS", "ImpulseResponse"};

figure
plot(irDataR.Time, irDataR.Amplitude)
title('R Measured')
grid on
xlim([0.1 0.2])
ylim([-0.04 0.04])

%%

binauralIrL = irDataL.Amplitude;
binauralIrR = irDataR.Amplitude;

% Compute spectrogram

% Generated by MATLAB(R) 9.13 and Signal Processing Toolbox 9.1.
% Generated on: 21-Mar-2023 15:15:16

% Parameters
frequencyLimits = [0 1]*pi; % Normalized frequency (rad/sample)
overlapPercent = 50;

%%
% Index into signal time region of interest
timeLimits = [1 176400]; % samples
binauralIrL_ROI = binauralIrL(:);
binauralIrL_ROI = binauralIrL_ROI(timeLimits(1):timeLimits(2));

limits = [-100 -50];

close all
% Compute spectral estimate
% Run the function call below without output arguments to plot the results
figure
pspectrum(binauralIrL_ROI, ...
    'spectrogram', ...
    'FrequencyLimits',frequencyLimits, ...
    'OverlapPercent',overlapPercent);
clim(limits)


% Compute spectrogram
% Index into signal time region of interest
timeLimits = [1 176400]; % samples
binauralIrR_ROI = binauralIrR(:);
binauralIrR_ROI = binauralIrR_ROI(timeLimits(1):timeLimits(2));

% Compute spectral estimate
% Run the function call below without output arguments to plot the results
figure
pspectrum(binauralIrR_ROI, ...
    'spectrogram', ...
    'FrequencyLimits',frequencyLimits, ...
    'OverlapPercent',overlapPercent);
clim(limits)

[RT,DRR,CTE,CFS,EDT] = irStats(['audio' filesep 'impulse_binaural.wav'], 'graph', true);

%%
audiowrite(['audio' filesep 'impulse2_binaural.wav'], [irDataL.Amplitude irDataR.Amplitude], fs)

audioFiles = {'musicDoubleBass', 'musicGuitar', 'musicSaxophone', 'musicViolin'};
numFiles = length(audioFiles);
for i = 1:numFiles
    audioFile = audioFiles{i};
    [audioStore, fs] = audioread(['sourceAudio/' audioFile '.wav']);
    if i == 1
        audio = audioStore;
    else
        audio = audio + audioStore;
    end
end
audioFile = 'music';
audioOutL = conv(audio, irDataL.Amplitude);
audioOutR = conv(audio, irDataR.Amplitude);
    
audiowrite(['audio' filesep audioFile '_binaural.wav'], [audioOutL audioOutR], fs)

audioFiles = {'piano', 'musicDoubleBass', 'musicGuitar', 'musicSaxophone', 'musicViolin'};
numFiles = length(audioFiles);

for i = 1:numFiles
    audioFile = audioFiles{i};
    [audio, fs] = audioread(['sourceAudio/' audioFile '.wav']);
    
    audioOutL = conv(audio, irDataL.Amplitude);
    audioOutR = conv(audio, irDataR.Amplitude);
    
    audiowrite(['audio' filesep audioFile '2_binaural.wav'], [audioOutL audioOutR], fs)
end

%%

clear all
close all

fs = 44.1e3;
ess = sweeptone(6, 4, fs);

load(['RIR' filesep 'PhDRoom.mat'])
irData24 = measurementData{"omni_MLS_24dBFS_C", "ImpulseResponse"};
irData6 = measurementData{"omni_MLS_6dBFS_C", "ImpulseResponse"};

load(['RIR' filesep 'BCCorridor_2ndOrder.mat'])
irData = measurementData{"omni_ESS", "ImpulseResponse"};

ir2 = irData.Amplitude;

load(['RIR' filesep 'BCCorridor_1stOrder.mat'])
irData = measurementData{"omni_ESS", "ImpulseResponse"};

ir1 = irData.Amplitude;

figure
plot(irData24.Time, irData24.Amplitude)
title('Omni24 Measured')

figure
plot(irData6.Time, irData6.Amplitude)
title('Omni6 Measured')

pathLength = 1;
c = 344;


sampleDelay = pathLength * fs / c;
fracDelay = sampleDelay - floor(sampleDelay);
sampleDelay = floor(sampleDelay);

ir = [zeros(sampleDelay, 1); (1 - fracDelay); fracDelay; 0];
tvec = (1:length(ir)) / fs;

figure
plot(tvec, ir)
title('Ref ir')
xlim([0 4])

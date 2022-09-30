% Create plots of how BTM changes with variables.

%% Create data

disp('Create data');

wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;
wedge = 190:10:350;
bendingAngle = 190:10:350;
minAngle = 0:10:180;

[result, geometry] = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, wedge, bendingAngle, minAngle);

result = rmfield(result,'i');

numResults = length(result);

%% Process data
disp('Process data');

wedge = [geometry.wedge];
bendingAngle = [geometry.bendingAngle];
minAngle = [geometry.minAngle];

tfmag = ([result.tfmag])';
fvec = result.fvec;

index = find(fvec(1,:) > 20000, 1);

dc = tfmag(:,1);
nyq = tfmag(:,index);
b0 = nyq - dc;

fcmag = dc - 3;
fcIndex = zeros(numResults, 1);

for i = 1:numResults
    fcIndex(i) = find(tfmag(i,:) < fcmag(i), 1);
end

fc = fvec(fcIndex)';

numOctaves = log2(20000 ./ fc);
gradient = (nyq + 3) ./ numOctaves;

data = [wedge, bendingAngle, minAngle, dc, nyq, b0, fc, gradient];

%% Save data
disp('Save data');

writematrix(data, 'BTMData');

%% Plots

for i = 1:numResults
    f = figure(1);
    movegui(f,'northeast');
    semilogx(result(i).fvec, result(i).tfmag)
    xlim([20, 20000])
    ylim([-40 0])
    text(fc(i), dc(i)-3, ['\leftarrow f_c: ', num2str(fc(i)), ' Hz'])
    text(fc(i), dc(i)-6, [num2str(gradient(i)), ' dB/octave'])
    title(['w: ', num2str(wedge(i)), ' bA: ', num2str(bendingAngle(i)), ' mA: ', num2str(minAngle(i))])
end


%% Data plots
close all

compare = [wedge, bendingAngle, minAngle, [geometry.source], wedge - [geometry.receiver]];

receiver = [geometry.receiver];

BTMDataPlot(receiver, wedge, bendingAngle, 'receiver', gradient, fc)

BTMDataPlot(bendingAngle, wedge, minAngle, 'Bending Angle', gradient, fc)
BTMDataPlot(minAngle, wedge, bendingAngle, 'Minimum Angle', gradient, fc)
BTMDataPlot(wedge, bendingAngle, minAngle, 'Wedge', gradient, fc)

toFind = [receiver, wedge];
[~, ~, idx] = unique(toFind, 'rows');

index = idx == 1;

rwfc = fc(index);
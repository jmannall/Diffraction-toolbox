%% Create result and geometry structs

close all

set(groot, 'defaultLineLineWidth', 1)



wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;

step = 10;
minw = 180;
maxw = 360;
shadowZone = true;

[result, geometry] = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, step, shadowZone, minw, maxw);

result = rmfield(result,'i');


%% Find unique indexes
input = transpose([[geometry.wedge]; [geometry.bendingAngle]; [geometry.minAngle]]);

[input, index, ~] = unique(input, 'rows');

result = result(index);
geometry = geometry(index);

bAindex = [geometry.bendingAngle] == 180;

result(bAindex) = [];
geometry(bAindex) = [];

%% Plot
numResults = length(result);

tfmag = zeros(numResults, 1);

for i = 1:numResults
    store = result(i).tfmag;
    tfmag(i) = store(44);
end

x = 1:numResults;
wedge = [geometry.wedge];
bendingAngle = [geometry.bendingAngle];
minAngle = [geometry.minAngle];
a = -8;
b = 25;
c = 800;
d = 57;
e = 1.5;
f = 1;
h = 1 / 9;
test = a * log10(b * wedge + c) .^ e + d - 10 * log10(h * bendingAngle) .^ 2;
test = a * log10(b * wedge + c) .^ e + d;
test = a * log10(b * wedge + c) .^ e + d - (wedge - 170) / 180 .*f .* log10(bendingAngle - 189) .^ 1.5 - 0.4 * log(minAngle + 1) .* (wedge - 189) / 360 .* (360 - bendingAngle) / 180;

tflin = 10 .^ (tfmag / 20);

close all

figure
plot(x,tflin)

figure
plot(x, tfmag)
hold on
plot(x, test)
hold off
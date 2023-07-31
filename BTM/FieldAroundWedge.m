function FieldAroundWedge(idx)

num = 500;
startIdx = idx * num;
endIdx = (idx + 1) * num;

fs = 96e3;
nfft = 8192;
c = 344;

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 0, 'noDirect', false);

wedgeLength = 20;
wedgeIndex = 270;
wedgeSize = 10;
x = cosd(wedgeIndex);
y = wedgeSize * sind(wedgeIndex);
corners = [0 0 0
    wedgeSize 0 0
    x y 0
    0 0 wedgeLength
    wedgeSize 0 wedgeLength
    x y wedgeLength
    wedgeSize y 0
    wedgeSize y wedgeLength];
        
planeCorners = [1 3 6 4
    3 7 8 6
    2 1 4 5
    4 6 8 5
    1 2 7 3
    2 5 8 7];

planeRigid = [1 0 1 0 0 0];

x = -3.05:0.1:3.05;
y = -3.05:0.1:3.05;
numPos = length(x) * length(y);

source = [-1, -2, wedgeSize / 2];
receiver = zeros(numPos, 3);

numFreq = nfft / 4;
template = -200 * ones(numFreq, numPos);
store = struct('complete', template, 'dir', template, 'spec', template, 'diff', template);

count = 0;
for i = x
    for j = y
        count = count + 1;
        receiver(count,:) = [i, j, wedgeSize / 2];
        if count > startIdx && count <= endIdx
            if i < 0 || j > 0
                [~, tfmag, ~, fvec] = SingleBTM(source, receiver(count,:), corners, planeCorners, planeRigid, controlparameters, false);
                store.complete(:,count) = tfmag.complete(1:numFreq);
                store.dir(:,count) = tfmag.direct(1:numFreq);
                store.spec(:,count) = tfmag.geom(1:numFreq);
                store.diff(:,count) = tfmag.diff1(1:numFreq);
            end
        end
    end
end

n = length(x);
xPlot = reshape(receiver(:,1), [], n);
yPlot = reshape(receiver(:,2), [], n);

template = -200 * ones(length(y), n, numFreq);
fields = struct('complete', template, 'dir', template, 'spec', template, 'diff', template);
for i = 1:numFreq
    fields.complete(:,:,i) = reshape(store.complete(i,:), [], n);
    fields.dir(:,:,i) = reshape(store.dir(i,:), [], n);
    fields.spec(:,:,i) = reshape(store.spec(i,:), [], n);
    fields.diff(:,:,i) = reshape(store.diff(i,:), [], n);
end
fvec = fvec(1:numFreq);

index = DataHash({source, x, y, wedgeIndex, wedgeLength, fs, nfft, c});
savePath = ['results' filesep 'fields'];
CheckFileDir(savePath);
fileName = [num2str(idx) '_' num2str(index)];
save([savePath filesep fileName], "fields", "xPlot", "yPlot", "fvec", '-v7.3')

% close all
% figure
% surf(xPlot, yPlot, fields.complete(:,:,1), 'FaceColor','interp', 'EdgeColor','none')
% view([0 90])
% clim([-50 0])
% colorbar
end

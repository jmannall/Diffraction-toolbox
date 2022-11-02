function [result, geometry] = SingleDiffractionGratingArray(numPillars, gratingWidth, pillarWidth, pillarHeight, receiver, controlparameters, createPlot)
    
    geometry = struct('numPillars', numPillars, 'gratingWidth', gratingWidth, 'pillarWidth', pillarWidth, 'pillarHeight', pillarHeight, 'receiver', receiver);
    pillarWidth = pillarWidth / 2;
    numReceivers = length(receiver);

    n = (numPillars - 1) / 2;
    pillarCenters = 3 * gratingWidth * [-fliplr(1:n), 0, 1:n];
    xDist = min(abs(receiver(:,2) - pillarCenters), [], 2) <= pillarWidth;
    yDist = abs(receiver(:,1)) <= pillarWidth;
    inPillar = xDist & yDist;

    toCalculate = numReceivers - sum(inPillar);
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({numPillars, gratingWidth, pillarWidth, pillarHeight, receiver, controlparameters});
    [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra, rtemplate] = BTMArrayFileHandling(mFile, index, numReceivers);
    template.complete = zeros(3, 1);
    template.direct = zeros(3, 1);
    template.geom = zeros(3, 1);
    template.diff1 = zeros(3, 1);
    template.diffHod = zeros(3, 1);

    saveCount = 1;    
    if resultExists
        load(loadPath, "result");
        % result = UniformResult(result);
        saveCount = CreateSaveCount(result, numReceivers, numSaves, filesPerSave);
    end
    if saveCount <= numSaves
        inPillari = ReshapeForParfor(inPillar, extra, filesPerSave);
        receiveri(:,1,:) = ReshapeForParfor(receiver(:,1), extra, filesPerSave);
        receiveri(:,2,:) = ReshapeForParfor(receiver(:,2), extra, filesPerSave);
        receiveri(:,3,:) = ReshapeForParfor(receiver(:,3), extra, filesPerSave);
        tic
        parfor j = saveCount:numSaves
            result = repmat(rtemplate, 1, 1);
            count = (j - 1) * filesPerSave;
            start = rem(count, filesPerSave);
            batchSize = min(filesPerSave, numReceivers - count);
            for i = (start + 1):batchSize
                inPillar = inPillari(i,j);
                receiver = receiveri(i,:,j);
                if inPillar
                    [tfmag, tfcomplex] = deal(repmat(template, 1, 1));

                    result(i).tfmag = tfmag
                    result(i).fvec = zeros(3, 1);
                    result(i).tfcomplex = tfcomplex;
                else
                    [result(i).tfmag, result(i).fvec, result(i).tfcomplex] = SingleDiffractionGrating(numPillars, gratingWidth, pillarWidth, pillarHeight, receiver, controlparameters, createPlot);
                end
                result(i).i = (j - 1) * filesPerSave + i;
            end
            savepathI = [savePath, '_', num2str(j)];
            ParforSaveArray(savepathI, result, geometry);
        end
        toc
        result = ProcessArrayResults(fileName, index, savePath, numSaves, geometry);
    end 
end
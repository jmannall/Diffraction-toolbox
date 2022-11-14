function [trainingData, targetData, fvec, fc, fidx, index, savePath] = CreateBtmTrainingData(numInputs, controlparameters, index)

    [geometry, trainingData] = RandomGeometryWedge(numInputs);

    % Create file info
    mFile = mfilename('fullpath');
    if nargin < 3
        index = DataHash({geometry, controlparameters});
    end
    [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra, rtemplate] = BTMArrayFileHandling(mFile, index, numInputs);

    saveCount = 1;
    if resultExists
        load(loadPath, "result");
        saveCount = CreateSaveCount(result, numInputs, numSaves, filesPerSave);
    end
    if saveCount <= numSaves
        wedgeIndexi = ReshapeForParfor(geometry.wedgeIndex, extra, filesPerSave);
        wedgeLengthi = ReshapeForParfor(geometry.wedgeLength, extra, filesPerSave);
        thetaSi = ReshapeForParfor(geometry.minAngle, extra, filesPerSave);
        thetaRi = ReshapeForParfor(geometry.minAngle + geometry.bendingAngle, extra, filesPerSave);
        rSi = ReshapeForParfor(geometry.rS, extra, filesPerSave);
        rRi = ReshapeForParfor(geometry.rR, extra, filesPerSave);
        zSi = ReshapeForParfor(geometry.truezS, extra, filesPerSave);
        zRi = ReshapeForParfor(geometry.truezR, extra, filesPerSave);
        tic
        parfor j = saveCount:numSaves
            controlparametersi = controlparameters;
            result = repmat(rtemplate, 1, 1);
            count = (j - 1) * filesPerSave;
            start = rem(count, filesPerSave);
            batchSize = min(filesPerSave, numInputs - count);
            for i = (start + 1):batchSize
                wedgeIndex = wedgeIndexi(i,j);
                wedgeLength = wedgeLengthi(i,j);
                thetaS = thetaSi(i,j);
                thetaR = thetaRi(i,j);
                rS = rSi(i,j);
                rR = rRi(i,j);
                zS = zSi(i,j);
                zR = zRi(i,j);

                L = sqrt((rS + rR) ^ 2 + (zR - zS) ^ 2);
                controlparametersi.Rstart = L;
                [result(i).tfmag, result(i).fvec, result(i).tfcomplex] = SingleWedgeInterpolated(wedgeLength,wedgeIndex,thetaS,thetaR,rS,rR,zS,zR,controlparametersi,false);
                result(i).i = (j - 1) * filesPerSave + i;
            end
            savepathI = [savePath, '_', num2str(j)];
            ParforSaveArray(savepathI, result, geometry);
        end
        toc
        result = ProcessArrayResults(fileName, index, savePath, numSaves, controlparameters, geometry);
    end
    tfmag = [result.tfmag];
    tfmag = max(min(tfmag, 128), -128);
    fvec = result(1).fvec;
    [targetData, fc, fidx] = CreateFrequencyNBands(tfmag, fvec, 12);
end
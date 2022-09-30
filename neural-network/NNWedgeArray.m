function [result, NNinput] = NNWedgeArray(geometry, NNinput, gParameters, numObservations, controlparameters)
    
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({gParameters, numObservations, controlparameters, geometry.distribution});
    [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra, rtemplate] = BTMArrayFileHandling(mFile, index, numObservations);

    saveCount = 1;    
    if resultExists
        load(loadPath, "result", "geometry", "NNinput");
        saveCount = CreateSaveCount(result, numInputs, numSaves);
    end
    if saveCount < numSaves
        wedgeLengthi = ReshapeForParfor(geometry.wedgeLength, extra, filesPerSave);
        wedgeIndexi = ReshapeForParfor(geometry.wedgeIndex, extra, filesPerSave);
        thetaSi = ReshapeForParfor(geometry.source, extra, filesPerSave);
        thetaRi = ReshapeForParfor(geometry.receiver, extra, filesPerSave);
        radiusSi = ReshapeForParfor(geometry.radiusS, extra, filesPerSave);
        radiusRi = ReshapeForParfor(geometry.radiusR, extra, filesPerSave);
        zSi = ReshapeForParfor(geometry.zS, extra, filesPerSave);
        zRi = ReshapeForParfor(geometry.zR, extra, filesPerSave);
        tic
        parfor j = saveCount:numSaves
            result = repmat(rtemplate, 1, 1);
            count = (j - 1) * filesPerSave;
            start = rem(count, filesPerSave);
            batchSize = min(filesPerSave, numInputs - count);
            for i = (start + 1):batchSize
                wedgeLength = wedgeLengthi(i,j);
                wedgeIndex = wedgeIndexi(i,j);
                thetaS = thetaSi(i,j);
                thetaR = thetaRi(i,j);
                radiusS = radiusSi(i,j);
                radiusR = radiusRi(i,j);
                zS = zSi(i,j);
                zR = zRi(i,j);

                [ir, tfmag, result(i).tvec, result(i).fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters, false);
                result(i).i = (j - 1) * filesPerSave + i;
                result(i).ir = ir.diff;
                result(i).tfmag = tfmag.diff;
                result(i).tfcomplex = tfcomplex.diff;
            end
            savepathI = [savePath, '_', num2str(j)];
            ParforSaveNNArray(savepathI, result, geometry, NNinput);
        end
        toc        
        result = ProcessArrayResults(fileName, index, savePath, numSaves, geometry);
    end 
end
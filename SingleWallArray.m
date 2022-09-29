function [result, geometry] = SingleWallArray(geometry, wallHeight, radiusS, radiusR, zS, zR, fs)
    
    numInputs = length(geometry.wallThickness);
    
    % Result template
    rtemplate.ir = [];
    rtemplate.tfmag = [];
    rtemplate.tvec = [];
    rtemplate.fvec = [];
    rtemplate.tfcomplex = [];
    rtemplate.i = [];
    
    result = repmat(rtemplate, 1, 1);
    
    index = DataHash({geometry, wallHeight, radiusS, radiusR, zS, zR, fs});
    
    % Create file info
    mFile = mfilename('fullpath');
    [inFilePath,fileStem] = fileparts(mFile);
    fileStem = [fileStem, '_', num2str(index)];
    savepath = ['results\', fileStem];
    loadpath = ['results\', fileStem];
    
    numPerSave = 50;
    numSaves = ceil(numInputs / numPerSave);
    saveIndex = 1:numSaves;

    i = 0;
    loadpathAll = [loadpath, '.mat'];
    test = exist(loadpathAll, "file");

    complete = false;
    if test == 0
        test = 2;
        while test == 2 && i < numSaves
            i = i + 1;
            loadpathI = [loadpath, '_', num2str(saveIndex(i)), '.mat'];
            test = exist(loadpathI, "file");
        end
        loadpath = [loadpath, '_', num2str(saveIndex(max(1,i - 1))), '.mat'];
    else
        loadpath = loadpathAll;
        complete = true;
        disp('Load from save');
    end

    test = exist(loadpath, "file");
    saveCount = 1;
    
    if test == 2
        load(loadpath, "result");
        count = result(end).i;
        saveCount = floor(count / numPerSave) + 1;
    end
    if ~complete
        extra = numPerSave - rem(numInputs, numPerSave);
        wallThickness = [geometry.wallThickness; zeros(extra, 1)];
        wallThicknessi = reshape(wallThickness,numPerSave,[]);
        thetaS = [geometry.source; zeros(extra, 1)];
        thetaSi = reshape(thetaS,numPerSave,[]);
        thetaR = [geometry.receiver; zeros(extra, 1)];
        thetaRi = reshape(thetaR,numPerSave,[]);
        tic
        parfor j = saveCount:numSaves
            result = repmat(rtemplate, 1, 1);
            count = (j - 1) * numPerSave;
            start = rem(count, numPerSave);
            batchSize = min(numPerSave, numInputs - count);
            for i = (start + 1):batchSize
                wallThickness = wallThicknessi(i,j);
                thetaS = thetaSi(i,j);
                thetaR = thetaRi(i,j);
                [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec, result(i).tfcomplex] = SingleWall(wallHeight,wallThickness,thetaS,thetaR,radiusS,radiusR,zS,zR,fs);
                result(i).i = i;
            end
            savepathI = [savepath, '_', num2str(saveIndex(j))];
            ParsaveSwa(savepathI, result, geometry);
        end
        toc
        resultAll = [];
        disp('Compiling saves. DO NOT STOP SCRIPT!')
        loadstem = ['results\', fileStem, '_'];
        for i = 1:numSaves
            loadpathI = [loadstem, num2str(saveIndex(i)), '.mat'];
            load(loadpathI, "result");
            resultAll = [resultAll, result];
        end
        result = resultAll;
        save(savepath, "result", "geometry", '-v7.3')
        for i = 1:numSaves
            loadpathI = [loadstem, num2str(saveIndex(i)), '.mat'];
            delete(loadpathI);
        end
    end 
end
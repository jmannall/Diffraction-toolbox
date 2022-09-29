function [result, NNinput] = NNWedgeArray(geometry, NNinput, gParameters, numObservations, fs)

    % Result template
    rtemplate.ir = [];
    rtemplate.tfmag = [];
    rtemplate.tvec = [];
    rtemplate.fvec = [];
    rtemplate.tfcomplex = [];
    rtemplate.i = [];
    
    index = DataHash({gParameters, numObservations, fs, geometry.distribution});
    
    % Create file info
    mFile = mfilename('fullpath');
    [inFilePath,fileStem] = fileparts(mFile);
    fileStem = [fileStem, '_', num2str(index)];
    savepath = ['results\', fileStem];
    loadpath = ['results\', fileStem];
    
    numPerSave = 50;
    numSaves = ceil(numObservations / numPerSave);
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
        load(loadpath, "result", "geometry", "NNinput");
        count = result(end).i;
        saveCount = floor(count / numPerSave) + 1;
    end
    if ~complete
        extra = numPerSave - rem(numObservations, numPerSave);
        wedgeLength = [geometry.wedgeLength; zeros(extra, 1)];
        wedgeLengthi = reshape(wedgeLength,numPerSave,[]);
        wedgeIndex = [geometry.wedgeIndex; zeros(extra, 1)];
        wedgeIndexi = reshape(wedgeIndex,numPerSave,[]);
        thetaS = [geometry.thetaS; zeros(extra, 1)];
        thetaSi = reshape(thetaS,numPerSave,[]);
        thetaR = [geometry.thetaR; zeros(extra, 1)];
        thetaRi = reshape(thetaR,numPerSave,[]);
        radiusS = [geometry.radiusS; zeros(extra, 1)];
        radiusSi = reshape(radiusS,numPerSave,[]);
        radiusR = [geometry.radiusR; zeros(extra, 1)];
        radiusRi = reshape(radiusR,numPerSave,[]);
        zS = [geometry.zS; zeros(extra, 1)];
        zSi = reshape(zS,numPerSave,[]);
        zR = [geometry.zR; zeros(extra, 1)];
        zRi = reshape(zR,numPerSave,[]);
        tic
        parfor j = saveCount:numSaves
            result = repmat(rtemplate, 1, 1);
            count = (j - 1) * numPerSave;
            start = rem(count, numPerSave);
            batchSize = min(numPerSave, numObservations - count);
            for i = (start + 1):batchSize
                count = (j - 1) * batchSize + i;
                wedgeLength = wedgeLengthi(i,j);
                wedgeIndex = wedgeIndexi(i,j);
                thetaS = thetaSi(i,j);
                thetaR = thetaRi(i,j);
                radiusS = radiusSi(i,j);
                radiusR = radiusRi(i,j);
                zS = zSi(i,j);
                zR = zRi(i,j);
%                 disp(['Wedge length: ', num2str(wedgeLength)]);
%                 disp(['WedgeIndex: ', num2str(wedgeIndex)]);
%                 disp(['thetaS: ', num2str(thetaS)]);
%                 disp(['thetaR: ', num2str(thetaR)]);
%                 disp(['zS: ', num2str(zS)]);
%                 disp(['zR: ', num2str(zR)]);
                [ir, tfmag, result(i).tvec, result(i).fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs);
                result(i).i = count;
                result(i).ir = ir.diff;
                result(i).tfmag = tfmag.diff;
                result(i).tfcomplex = tfcomplex.diff;
            end
            savepathI = [savepath, '_', num2str(saveIndex(j))];            
            ParsaveNNwa(savepathI, result, geometry, NNinput);
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
        save(savepath, "result", "geometry", "NNinput", '-v7.3')
        for i = 1:numSaves
            loadpathI = [loadstem, num2str(saveIndex(i)), '.mat'];
            delete(loadpathI);
        end
    end 
end
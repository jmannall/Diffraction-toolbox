function [trainingData, targetData, DC, fvec, fc, fidx, index, savePath, tfmag, tfmagI] = CreateBtmTrainingData(numInputs, controlparameters, index)

    geometry = RandomGeometryWedge_Run2(numInputs);

    % Create file info
    mFile = mfilename('fullpath');
    saveAll = false;
    if nargin < 3
        index = DataHash({geometry, controlparameters});
        disp(['DataHash: ' num2str(index)])
        saveAll = true;
    end
    [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra] = BTMArrayFileHandling(mFile, index, numInputs);

    fs = controlparameters.fs / 2;
    nfft = controlparameters.nfft / 2;
    fvec = fs/nfft*[0:nfft/2-1];

    NNSavePath = [savePath, '-', num2str(round(fs / 1e3)), '_NN'];
    NNSaveExists = exist([cd, filesep, NNSavePath, '.mat'], "file");

    if NNSaveExists == 2
        %load([NNSavePath, '.mat'], "trainingData", "targetData", "fvec", "fc", "fidx", "index", "savePath");
        load([NNSavePath, '.mat']);
        DC = 1;
        % disp('Load NN data')
    else
        rtemplate.tfmag = [];
        rtemplate.i = [];

        saveCount = 1;
        if resultExists
            load(loadPath, "result", "geometry");
            saveCount = CreateSaveCount(result, numInputs, numSaves, filesPerSave);
        end
        if saveCount <= numSaves
            wedgeIndexi = ReshapeForParfor(geometry.wedgeIndex, extra, filesPerSave);
            wedgeLengthi = ReshapeForParfor(geometry.wedgeLength, extra, filesPerSave);
            thetaSi = ReshapeForParfor(geometry.minAngle, extra, filesPerSave);
            thetaRi = ReshapeForParfor(geometry.minAngle + geometry.bendingAngle, extra, filesPerSave);
            rSi = ReshapeForParfor(geometry.rS, extra, filesPerSave);
            rRi = ReshapeForParfor(geometry.rR, extra, filesPerSave);
            zSi = ReshapeForParfor(geometry.zS, extra, filesPerSave);
            zRi = ReshapeForParfor(geometry.zR, extra, filesPerSave);
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

                    result(i).tfmagI = SingleWedgeInterpolated(wedgeLength,wedgeIndex,thetaS,thetaR,rS,rR,zS,zR,controlparametersi,false);
                    if saveAll
                        [~, tempVar] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,rS,rR,zS,zR,controlparametersi,false);
                        result(i).tfmag = tempVar.diff1;
                    end
                    result(i).i = i;
                end
                savepathI = [savePath, '_', num2str(j)];
                ParforSaveArray(savepathI, result, geometry);
            end
            toc
            result = ProcessArrayResults(fileName, index, savePath, numSaves, controlparameters, geometry);
        end
        tfmagI = [result.tfmagI];
        tfmagI = tfmagI(1:end / 2,:);
        tfmagI = max(min(tfmagI, 128), -128);
        DC = tfmagI(1,:);
        [targetData, fc, fidx] = CreateFrequencyNBands(tfmagI, fvec, 8);
        trainingData = CreateNNinput(geometry);

        if saveAll
            tfmag = [result.tfmag];
            tfmag = tfmag(1:end / 2,:);
            tfmag = max(min(tfmag, 128), -128);
            save(NNSavePath, "trainingData", "targetData", "DC", "tfmagI", "tfmag", "fvec", "fc", "fidx", "index", "savePath", '-v7.3')
            save([NNSavePath, '_small'], "trainingData", "targetData", "DC", "fvec", "fc", "fidx", "index", "savePath", '-v7.3')
        else
            save(NNSavePath, "trainingData", "targetData", "DC", "fvec", "fc", "fidx", "index", "savePath", '-v7.3')
        end
        delete([savePath, '.mat'])
    end
end
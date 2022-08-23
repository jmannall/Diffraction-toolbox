function [result, input] = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, step, shadowZone, minw, maxw)

    input = Geometry(step, shadowZone, minw, maxw);
    
    if isa(input, 'double')
        result = 0;
        return
    end
    
    numInputs = length(input);
    
    % Result template
    rtemplate.ir = [];
    rtemplate.tfmag = [];
    rtemplate.tvec = [];
    rtemplate.fvec = [];
    rtemplate.tfcomplex = [];
    rtemplate.i = [];
    
    result = repmat(rtemplate, 1, 1);
    
    index1 = DataHash(input);
    index2 = DataHash([wedgeLength, radiusS, radiusR, zS, zR, fs]);
    
    index = [index1, index2];
    
    % Create file info
    mFile = mfilename('fullpath');
    [inFilePath,fileStem] = fileparts(mFile);
    fileStem = [fileStem, '_', num2str(index)];
    savepath = ['results\', fileStem];
    loadpath = ['results\', fileStem, '.mat'];
    
    test = exist(loadpath, "file");
    count = 0;
    complete = false;
    
    if test == 2
        load(loadpath, "result");
        count = result(end).i;
        if count == numInputs
            complete = true;
        end
        save(savepath, "result");
    end
    if ~complete
        for i = (count + 1):numInputs
            wedgeIndex = input(i).wedge;
            thetaS = input(i).source + 0.01;
            thetaR = input(i).receiver - 0.01;
            if thetaR == 359.99
                thetaR = 359.98;
            end
            if wedgeIndex == 180
                wedgeIndex = 180.01;
            elseif wedgeIndex == 360
                wedgeIndex = 359.99;
            end
            [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec, result(i).tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs);
            result(i).i = i;
            save(savepath, "result");
        end
    end
end
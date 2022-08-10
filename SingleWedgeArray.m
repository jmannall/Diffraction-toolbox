function result = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, step, shadowZone)

input = Geometry(step, shadowZone);

numInputs = length(input);

% Result template
rtemplate.ir = [];
rtemplate.tfmag = [];
rtemplate.tvec = [];
rtemplate.fvec = [];
rtemplate.tfcomplex = [];

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

if test == 2
    load(loadpath, "result");
else
    for i = 1:numInputs
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
    end
end

save(savepath, "result");

end
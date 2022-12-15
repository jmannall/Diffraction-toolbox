function tfmag = CreateReferenceNNData(inputData, controlparameters)
    numData = size(inputData, 2);

    thetaS = rad2deg(inputData(3,:));
    thetaR = rad2deg(inputData(2,:) + inputData(3,:));
    wedgeIndex = rad2deg(inputData(1,:));
    rS = inputData(5,:);
    rR = inputData(6,:);
    zS = inputData(7,:);
    zR = inputData(8,:);
    wedgeLength = inputData(4,:);
    for i = 1:numData
        [tfmag, fvec] = SingleWedgeInterpolated(wedgeLength(i), wedgeIndex(i), thetaS(i), thetaR(i), rS(i), rR(i), zS(i), zR(i), controlparameters, false);
    end
    
    tfmag = CreateFrequencyNBands(tfmag, fvec, 12);
end
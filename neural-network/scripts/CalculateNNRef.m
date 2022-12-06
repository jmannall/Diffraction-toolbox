function [tfcomplex, ir] = CalculateNNRef(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, validPath, pathLengthDir, pathLengthNN)

    zS = zS + wedgeLength;
    zR = zR + wedgeLength;
    wedgeLength = 2 * wedgeLength;

    numReceivers = length(thetaR);
    tfcomplexRef = zeros(controlparameters.nfft / 2, numReceivers);
    for i = 1:numReceivers
        [~, ~, tfcomplexRef(:,i), ir{i}] = SingleWedgeInterpolated(wedgeLength,wedgeIndex,thetaS,thetaR(i),rS,rR(i),zS,zR(i),controlparameters,false);
    end
    
    tfcomplex(:, validPath) = tfcomplexRef(:, validPath) ./ pathLengthDir(validPath)';
    tfcomplex(:, ~validPath) = tfcomplexRef(:, ~validPath) ./ pathLengthNN(~validPath)';
end
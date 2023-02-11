function [tfcomplex, ir] = CalculateNNRef(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, dirVisible, pathLengthDir, pathLengthEd, doReflection)

    if doReflection
        zS = zS + wedgeLength;
        zR = zR + wedgeLength;
        wedgeLength = 2 * wedgeLength;
    end

    numReceivers = length(thetaR);
    tfcomplexRef = zeros(controlparameters.nfft / 2, numReceivers);
    [~, ~, tfcomplexRef, ir] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, false);

%     for i = 1:numReceivers
%         [~, ~, tfcomplexRef(:,i), ir] = SingleWedgeInterpolated(wedgeLength,wedgeIndex,thetaS,thetaR(i),rS,rR(i),zS,zR(i),controlparameters,false);
%     end
    
    tfcomplex(:, dirVisible) = tfcomplexRef(:, dirVisible) ./ pathLengthDir(dirVisible)';
    tfcomplex(:, ~dirVisible) = tfcomplexRef(:, ~dirVisible) ./ pathLengthEd(~dirVisible)';
end
function [tfmag, fvec, tfcomplex] = SingleUTDExtDaisyChain(data, phii, controlparameters)

    rS = [data.rS, data.W];
    rR = [data.rR, fliplr(data.W)];

    cumRs = cumsum(rS)';
    cumRr = fliplr(cumsum(rR))';

    wedgeIndex = data.wedgeIndex;
    numEdges = length(wedgeIndex);
    
    thetaS = [data.thetaS, zeros(1, numEdges - 1)];
    thetaR = [wedgeIndex(1:end - 1), data.thetaR];
    
    tfcomplex = zeros(4, numEdges + 1);

    if controlparameters.interpolated
        minAngle = min(thetaS, wedgeIndex - thetaR);
        bendingAngle = thetaR - thetaS;
        thetaS = minAngle;
        thetaR = minAngle + bendingAngle;
        UTDWedge = @(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii)SingleUTDWedgeInterpolated(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters);
    else
        UTDWedge = @(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii)SingleUTDWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters);
    end
    for i = 1:numEdges
        [~, fvec, tfcomplexStore] = UTDWedge(thetaS(i), thetaR(i), cumRs(i), cumRr(i), wedgeIndex(i), phii);
        tfcomplex(:,i) = (cumRs(i) + cumRr(i)) * tfcomplexStore;
    end
    

    tfcomplex(:,end) = 0.5^(numEdges - 1) .* (1 / data.L) .* prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
function [tfmag, fvec, tfcomplex] = SingleUTDApexDaisyChainV2(data, phi, always, c, controlparameters, withCorrection)

    numEdges = size(data.wedgeIndex, 1);

    L = data.L;
    rS = data.rS;
    rR = data.rR;
    W = data.W;

    rS = [rS, W];
    rR = [W, rR];
    wedgeIndex = data.wedgeIndex';
    
    epsilon = 1e-5;
    thetaS = [data.thetaS, epsilon * ones(1, numEdges)];
    thetaR = [wedgeIndex - epsilon, data.thetaR];

    scale = KimCorrection(data, numEdges, withCorrection);
    
    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);
    controlparameters.difforder = 1;
    for i = 1:numEdges
        [~, fvec, tfcomplexStore] = SingleUTDWedge(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phi, always, controlparameters, c);
        tfcomplex(:,i) = scale(i) * tfcomplexStore;
    end
    tfcomplex(:,end) = 0.5^(numEdges - 1) * prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
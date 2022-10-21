function [tfmag, fvec, tfcomplex] = UTDSingleApexDaisyChain(data, phii, controlparameters, withCorrection)

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

    [scale, A] = KimCorrection(data, numEdges, withCorrection);
    
    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);
    controlparameters.difforder = 1;
    for i = 1:numEdges
        [~, fvec, tfcomplexStore] = UTDSingleWedge(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phii, controlparameters);
        tfcomplex(:,i) = (rS(i) + rR(i)) * scale(i) * tfcomplexStore;
    end
    tfcomplex(:,end) = 0.5^(numEdges - 1) * (1 / L) * prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
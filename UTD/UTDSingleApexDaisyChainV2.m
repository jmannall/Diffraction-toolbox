function [tfmag, fvec, tfcomplex] = UTDSingleApexDaisyChainV2(data, phii, controlparameters, withCorrection)

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

    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);
    controlparameters.difforder = 1;

    if withCorrection
        [A, B] = KimCorrection(data, numEdges);    
        for i = 1:numEdges
            [~, fvec, tfcomplexStore] = UTDSingleWedge(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phii, controlparameters, A(i), B(i));
            tfcomplex(:,i) = (1 / sqrt(A(i) * B(i))) * tfcomplexStore;
        end
        c = controlparameters.c;
        k = 2 * pi * fvec/ c;
        front = exp(-1i .* k .* L)';
    else
        for i = 1:numEdges
            [~, fvec, tfcomplexStore] = UTDSingleWedge(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phii, controlparameters);
            tfcomplex(:,i) = scale(i) * tfcomplexStore;
        end
        front = 1;
    end

    tfcomplex(:,end) = front .* 0.5^(numEdges - 1) * prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
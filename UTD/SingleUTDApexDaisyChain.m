function [tfmag, fvec, tfcomplex] = SingleUTDApexDaisyChain(data, phii, controlparameters, withCorrection)

    rS = [data.rS, data.W];
    rR = [data.W, data.rR];
    wedgeIndex = data.wedgeIndex;
    numEdges = length(wedgeIndex);
    
    thetaS = [data.thetaS, zeros(1, numEdges)];
    thetaR = [wedgeIndex, data.thetaR];
    
    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);
    if withCorrection
        [A, B] = KimCorrection(data, numEdges);    
        for i = 1:numEdges
            [~, fvec, tfcomplexStore] = SingleUTDWedge(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phii, controlparameters, A(i), B(i));
            tfcomplex(:,i) = (1 / sqrt(A(i) * B(i))) * tfcomplexStore;
        end
        c = controlparameters.c;
        k = 2 * pi * fvec/ c;
        front = exp(-1i .* k .* data.L);
    else
        for i = 1:numEdges
            [~, fvec, tfcomplexStore] = SingleUTDWedge(thetaS(i), thetaR(i), rS(i), rR(i), wedgeIndex(i), phii, controlparameters);
            tfcomplex(:,i) = (rS(i) + rR(i)) * tfcomplexStore;
        end
        front = 1;
    end

    tfcomplex(:,end) = front' .* 0.5^(numEdges - 1) .* (1 / data.L) .* prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
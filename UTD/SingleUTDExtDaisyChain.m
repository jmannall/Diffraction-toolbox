function [tfmag, fvec, tfcomplex] = SingleUTDExtDaisyChain(data, phii, controlparameters)

    rS = [data.rS, data.W];
    rR = [data.rR, fliplr(data.W)];

    cumRs = cumsum(rS)';
    cumRr = fliplr(cumsum(rR))';

    wedgeIndex = data.wedgeIndex;
    numEdges = length(wedgeIndex);
    
    thetaS = [data.thetaS, zeros(1, numEdges)];
    thetaR = [wedgeIndex, data.thetaR];
    
    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);

    for i = 1:numEdges
        [~, fvec, tfcomplexStore] = SingleUTDWedge(thetaS(i), thetaR(i), cumRs(i), cumRr(i), wedgeIndex(i), phii, controlparameters);
        tfcomplex(:,i) = (cumRs(i) + cumRr(i)) * tfcomplexStore;
    end
    

    tfcomplex(:,end) = 0.5^(numEdges - 1) .* (1 / data.L) .* prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
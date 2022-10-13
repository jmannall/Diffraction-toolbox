function [ir, tfmag, tvec, fvec, tfcomplex] = SingleNthOrderBarrier1storder(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, numEdges, createPlot)

    [corners, planeCorners, ~, source, receiver] = CreateNthOrderBarrierGeometry(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, numEdges);
    [thetaS, thetaR, A, B, L, radiusS, radiusR, wedgeIndex, source, receiver] = CreateNthOrderBarrierGeometry1stOrder(corners, source, receiver, numEdges);
    
    if createPlot
        PlotGeometry(corners, planeCorners, source, receiver)
    end

    nir = 0;
    [irStore, tvecStore] = deal(cell(1, numEdges));
    [tfmag, tfcomplex] = deal(zeros(controlparameters.nfft / 2, numEdges + 1));

    radiusR = 10^6;
    radiusS = 1;
    controlparameters.Rstart = radiusR;
    for i = 1:numEdges
%         [irAll, tfmagAll, tvecStore{i}, fvec, tfcomplexAll] = SingleWedge(barrierHeight, wedgeIndex(i), thetaS(i), thetaR(i), radiusS(i), radiusR(i), zS, zR, controlparameters, createPlot);
        [irAll, tfmagAll, tvecStore{i}, fvec, tfcomplexAll] = SingleWedge(barrierHeight, wedgeIndex(i), thetaS(i), thetaR(i), radiusR, radiusS, zS, zR, controlparameters, createPlot);

        nir = max(nir, length(irAll.diff1));
        irStore{i} = radiusR * irAll.diff1;
        %tfmag(:,i) = tfmagAll.diff1;
        scale = (1 / sqrt(A(i) * B(i)));
        tfcomplex(:,i) = radiusR * scale * tfcomplexAll.diff1;
%         tfcomplex(:,i) = (1 / sqrt(A(i) * B(i))) * tfcomplexAll.diff1;
    end

    % Make all the impulse responses the same length
    ir = zeros(nir, numEdges);
    tvec = zeros(1, nir);
    for i = 1:numEdges
        idx = 1:length(irStore{i});
        ir(idx,i) = irStore{i};
        tvec(1,idx) = tvecStore{i};
    end
    
    % 0.5^M per connected wegdes 
    tfcomplex(:,end) = 0.5^(numEdges - 1) * (1 / L) * prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
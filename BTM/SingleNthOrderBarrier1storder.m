function [ir, tfmag, tvec, fvec, tfcomplex] = SingleNthOrderBarrier1storder(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, numEdges, createPlot)

    [corners, planeCorners, ~, source, receiver] = CreateNthOrderBarrierGeometry(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, numEdges);
    [thetaS, thetaR, radiusS, radiusR, wedgeIndex, source, receiver] = CreateNthOrderBarrierGeometry1stOrder(corners, source, receiver, numEdges);
    
    if createPlot
        PlotGeometry(corners, planeCorners, source, receiver)
    end

    nir = 0;
    [irStore, tvecStore] = deal(cell(1, numEdges));
    [tfmag, tfcomplex] = deal(zeros(controlparameters.nfft / 2, numEdges + 1));
    for i = 1:numEdges
        [irAll, tfmagAll, tvecStore{i}, fvec, tfcomplexAll] = SingleWedge(barrierHeight, wedgeIndex(i), thetaS(i), thetaR(i), radiusS(i), radiusR(i), zS, zR, controlparameters, createPlot);
        
        nir = max(nir, length(irAll.diff1));
        irStore{i} = irAll.diff1;
        tfmag(:,i) = tfmagAll.diff1;
        tfcomplex(:,i) = tfcomplexAll.diff1;
    end

    % Make all the impulse responses the same length
    ir = zeros(nir, numEdges);
    tvec = zeros(1, nir);
    for i = 1:numEdges
        idx = 1:length(irStore{i});
        ir(idx,i) = irStore{i};
        tvec(1,idx) = tvecStore{i};
    end
    
    tfcomplex(:,end) = prod(tfcomplex(:,1:numEdges), 2);
    tfmag(:,end) = mag2db(abs(tfcomplex(:,end)));
end
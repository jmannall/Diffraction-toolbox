function [tfmag, fvec, tfcomplex] = SingleBTMApexDaisyChainV3(source, receiver, apex, corners, planeCorners, controlparameters, data, withCorrection, createPlot)

    numEdges = size(apex, 1);
    % Create virtual sources and receivers
    vSource = [source; apex(1:numEdges - 1,:)];
    vReceiver = [apex(2:numEdges,:); receiver];

    L = data.L;

    rS = data.rS;
    rR = data.rR;
    W = data.W;

    rS = [rS, W];
    rR = [W, rR];

    scale = KimCorrection(data, numEdges, withCorrection);
    
    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);
    controlparameters.difforder = 1;
    for i = 1:numEdges
        vPlaneRigid = [zeros(1, i - 1), ones(1, 2), zeros(1, numEdges + 3 - i)];
        [~, ~, ~, fvec, tfcomplexStore] = SingleBTM(vSource(i,:), vReceiver(i,:), corners, planeCorners, vPlaneRigid, controlparameters, createPlot);
        tfcomplex(:,i) = scale(i) * tfcomplexStore.diff1;
    end
    tfcomplex(:,end) = 0.5^(numEdges - 1) * prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
function [tfmag, fvec, tfcomplex] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, rS, rR, L)

    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);
    controlparameters.difforder = 1;
    for i = 1:numEdges
        vPlaneRigid = [zeros(1, i - 1), ones(1, 2), zeros(1, numEdges + 3 - i)];
        [~, ~, ~, fvec, tfcomplexStore] = SingleBTM(vSource(i,:), vReceiver(i,:), corners, planeCorners, vPlaneRigid, controlparameters, false);
        tfcomplex(:,i) = (rS(i) + rR(i)) * tfcomplexStore.diff1;
    end

    tfcomplex(:,end) = 0.5^(numEdges - 1) * (1 / L) * prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
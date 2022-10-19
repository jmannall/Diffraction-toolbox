function [tfmag, fvec, tfcomplex] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, normReceiver, withCorrection, createPlot)

    numEdges = size(apex, 1);
    % Create virtual sources and receivers
    vSource = [source; apex(1:numEdges - 1,:)];
    vReceiver = [apex(2:numEdges,:); receiver];

    const = 10^6;
    vCorners = [corners(2:numEdges + 1,1:2), vSource(:,3)];
    vector = vSource - vCorners;
    vector = const * vector ./ vecnorm(vector, 2, 2);
    vSource = vCorners + vector;
    controlparameters.Rstart = const;

    if normReceiver
        vector = vReceiver - vCorners;
        vector = vector ./ vecnorm(vector, 2, 2);
        vReceiver = vCorners + vector;
    end

    scale = KimCorrection(data, numEdges, withCorrection);

    tfcomplex = zeros(controlparameters.nfft / 2, numEdges + 1);
    controlparameters.difforder = 1;
%     vSource(1,:) = source;
    for i = 1:numEdges
%         if i == 1
%             controlparameters.Rstart = data.rS;
%         else
%             controlparameters.Rstart = const;
%         end
        vPlaneRigid = [zeros(1, i - 1), ones(1, 2), zeros(1, numEdges + 3 - i)];
        [~, ~, ~, fvec, tfcomplexStore] = SingleBTM(vSource(i,:), vReceiver(i,:), corners, planeCorners, vPlaneRigid, controlparameters, createPlot);
        tfcomplex(:,i) = scale(i) * tfcomplexStore.diff1;
    end
    tfcomplex(:,end) = 0.5^(numEdges - 1) * (1 / data.L) * prod(tfcomplex(:,1:numEdges), 2);
    tfmag = mag2db(abs(tfcomplex));
end
function [tfmag, fvec, tfcomplex] = SingleBTMExtDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data)

    numEdges = size(apex, 1);

    % Create virtual sources and receivers
    vSource = [source; apex(1:numEdges - 1,:)];
    vReceiver = [apex(2:numEdges,:); receiver];

    rS = [data.rS, data.W];
    rR = [data.rR, fliplr(data.W)];

    cumRs = cumsum(rS)';
    cumRr = fliplr(cumsum(rR))';

    vCorners = [corners(2:numEdges + 1,1:2), vSource(:,3)];
    vector = vSource - vCorners;
    vector = cumRs .* vector ./ vecnorm(vector, 2, 2);
    vSource = vCorners + vector;

    vector = vReceiver - vCorners;
    vector = cumRr .* vector ./ vecnorm(vector, 2, 2);
    vReceiver = vCorners + vector;

    [tfmag, fvec, tfcomplex] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, cumRs, cumRr, data.L);
end

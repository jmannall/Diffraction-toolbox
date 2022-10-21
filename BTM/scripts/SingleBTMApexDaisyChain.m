function [tfmag, fvec, tfcomplex] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data)

    numEdges = size(apex, 1);

    % Create virtual sources and receivers
    vSource = [source; apex(1:numEdges - 1,:)];
    vReceiver = [apex(2:numEdges,:); receiver];

    rS = [data.rS, data.W];
    rR = [data.W, data.rR];
    
    [tfmag, fvec, tfcomplex] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, rS, rR, data.L);
end
function [tfmag, fvec, tfcomplex] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data)

    numEdges = size(apex, 1);

    % Create virtual sources and receivers
    vSource = [source; apex(1:numEdges - 1,:)];
    vReceiver = [apex(2:numEdges,:); receiver];

    rS = zeros(numEdges, 1);
    rR = ones(numEdges, 1);
    
    const = 1e7;
    vCorners = [corners(2:numEdges + 1,1:2), vSource(:,3)];
    vector = vSource - vCorners;
    vector = const .* vector ./ vecnorm(vector, 2, 2);
    vSource = vCorners + vector;
    controlparameters.Rstart = const;

    % In some geometries needed so that receivers can see the scattering
    % geometry
    if data.W < 1
        vCorners = [corners(2:numEdges + 1,1:2), vSource(:,3)];
    else
        vCorners = [apex(:,1:2), vSource(:,3)];
    end
    vector = vReceiver - vCorners;
    vector = rR .* vector ./ vecnorm(vector, 2, 2);
    vReceiver = vCorners + vector;

    [tfmag, fvec, tfcomplex] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, rS, (rR + const - 1) / const, data.L);
end

function [corners, planeCorners, source, receiver] = CreateDiffractionGratingGeometry(numPillars, gratingWidth, pillarHeight, receiver)

    if nargin < 4
        receiver = [0 5 5];
    end
    receiver(1) = abs(receiver(1));
    [corners, planeCorners, planeRigid] = CreatePillar(1, 1, gratingWidth, pillarHeight);
    for i = 2:(numPillars + 1) / 2
        for j = 0:1
            [cornersStore, planeCornersStore, planeRigidStore] = CreatePillar(i, j, gratingWidth, pillarHeight);
            corners = [corners; cornersStore];
            planeCorners = [planeCorners; planeCornersStore];
            planeRigid = [planeRigid, planeRigidStore];
        end
    end
    
    source = [0 -7 5];
end

    
    
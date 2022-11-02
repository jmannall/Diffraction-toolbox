function [corners, planeCorners, source, receiver] = CreateDiffractionGratingGeometry(numPillars, gratingWidth, pillarWidth, pillarHeight, receiver)

    if nargin < 5
        receiver = [0 5 5];
    end
    [corners, planeCorners, planeRigid] = CreatePillar(1, 1, gratingWidth, pillarWidth, pillarHeight);
    for i = 2:(numPillars + 1) / 2
        for j = 0:1
            [cornersStore, planeCornersStore, planeRigidStore] = CreatePillar(i, j, gratingWidth, pillarWidth, pillarHeight);
            corners = [corners; cornersStore];
            planeCorners = [planeCorners; planeCornersStore];
            planeRigid = [planeRigid, planeRigidStore];
        end
    end
    
    source = [-4.4 / 3.5 0 pillarHeight / 2];
end

    
    
function [corners, planeCorners, planeRigid, source, receiver] = CreateNthOrderBarrierGeometry(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, numDiffEdges)
    
    numEdges = numDiffEdges + 2;
    angles = linspace(0, 180, numDiffEdges)';
    
    corners = [-2 * radiusS -barrierRadius 0
        barrierRadius * sind(angles)  -barrierRadius * cosd(angles) zeros(numDiffEdges, 1)
        -2 * radiusR barrierRadius 0
        -2 * radiusS -barrierRadius barrierHeight
        barrierRadius * sind(angles) -barrierRadius * cosd(angles) barrierHeight * ones(numDiffEdges, 1)
        -2 * radiusR barrierRadius barrierHeight];
    
    planeCorners = zeros(numEdges + 2, numEdges);
    for i = 1:numEdges - 1
        planeCorners(i,1:4) = [i i + 1 i + 1 + numEdges i + numEdges];
    end
    planeCorners(numEdges, 1:4) = [numEdges, 1, numEdges + 1 2 * numEdges];
    planeCorners(numEdges + 1,:) = [1, fliplr(2:numEdges)];
    planeCorners(numEdges + 2,:) = numEdges + 1:2 * numEdges;
    
    planeRigid = [ones(1, numDiffEdges + 1), 0, 0, 0, 0];
    
    source = [-radiusS * cosd(thetaS), -radiusS * sind(thetaS) - barrierRadius, zS];
    receiver = [radiusR * sind(thetaR), -radiusR * cosd(thetaR) + barrierRadius, zR];
end
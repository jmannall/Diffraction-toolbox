function [source, receiver, Q, apex, corners, planeCorners, planeRigid, data] = GenerateNthOrderPath(numEdges, height)
    
    valid = false;
    while ~valid
        
        shift = 3;
        scale = 10 - shift;
        rS = shift + scale * rand(1);
        rR = shift + scale * rand(1);
        W = shift + scale * rand(1, numEdges - 1);
    
        wedgeIndex = 185 + (355 - 185) * rand(1, numEdges);
        thetaS = 0.1 + (wedgeIndex(1) - 180) * rand(1);
        thetaR = 180 + (wedgeIndex(end) - 180.1) * rand(1);

        data = struct('rS', rS, 'rR', rR, 'W', W, 'L', rS + sum(W) + rR, 'thetaS', thetaS, 'thetaR', thetaR, 'wedgeIndex', wedgeIndex);
    
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, rS, rR, W, height);
    end
end
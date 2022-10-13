function [thetaS, thetaR, radiusS, radiusR, wedgeIndex, source, receiver] = CreateNthOrderBarrierGeometry1stOrder(corners, source, receiver, numEdges)
    
    zS = source(:, 3);
    zR = receiver(:, 3);
    corners = corners(1:end / 2, [1,2]);
    source = source(:, [1,2]);
    receiver = receiver(:, [1,2]);
    
    wedgeIndex = zeros(1, numEdges);
    for i = 1:numEdges
        a = corners(i + 1, :) - corners(i, :);
        b = corners(i + 1, :) - corners(i + 2, :);
        wedgeIndex(i) = 360 - AngleBetweenVectors(a, b);
    end
    
    [virtualReceiver, virtualSource] = deal(zeros(numEdges - 1, 2));
    for i = 2:numEdges
        virtualSource(i - 1, :) =  corners(i, :);
        virtualReceiver(i - 1, :) =  corners(i + 1, :);
    end

    source = [source; virtualSource];
    receiver = [virtualReceiver; receiver];
    
    radiusS = vecnorm((corners(2:end - 1,:) - source)');
    radiusR = vecnorm((corners(2:end - 1,:) - receiver)');

%     radiusS = [vecnorm((corners(2,:) - source), 2), max(0.5, min(vecnorm((corners(3:end - 1,:) - virtualSource)'), 2))];
%     radiusR = [max(0.5, min(vecnorm((corners(2:end - 2,:) - virtualReceiver)'), 2)), vecnorm((corners(end - 1,:) - receiver), 2)];
%     
%     source = [source; virtualSource];
%     receiver = [virtualReceiver; receiver];

    thetaS = max(AngleBetweenVectors((corners(2:end - 1,:) - source)', (corners(2:end - 1,:) - corners(1:end - 2,:))'), 0.001);
    thetaR = wedgeIndex - max(AngleBetweenVectors((corners(2:end - 1,:) - receiver)', (corners(2:end - 1,:) - corners(3:end,:))'), 0.001);

    source = [source, zS * ones(numEdges, 1)];
    receiver = [receiver, zR * ones(numEdges, 1)];
end
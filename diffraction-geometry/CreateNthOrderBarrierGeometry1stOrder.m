function [thetaS, thetaR, A, B, L, rS, rR, wedgeIndex, source, receiver] = CreateNthOrderBarrierGeometry1stOrder(corners, source, receiver, numEdges)
    
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

    % Distances between wedges in the path
    [W, p] = deal(zeros(1,numEdges - 1));
    for i = 1:numEdges - 1
        W(i) = vecnorm((corners(i + 1,:) - corners(i + 2,:))');
    end
    
    radiusS = vecnorm((corners(2,:) - source)');
    radiusR = vecnorm((corners(end - 1,:) - receiver)');

    source = [source; virtualSource];
    receiver = [virtualReceiver; receiver];
    thetaS = max(AngleBetweenVectors((corners(2:end - 1,:) - source)', (corners(2:end - 1,:) - corners(1:end - 2,:))'), 0.001);
    thetaR = wedgeIndex - max(AngleBetweenVectors((corners(2:end - 1,:) - receiver)', (corners(2:end - 1,:) - corners(3:end,:))'), 0.001);


    % Total path length
    L = radiusS + sum(W) + radiusR;
    % Q is position of wedge
    % A = distance along path from Q to R * distance along path from Q to S / L
%     Atest = ([(sum(W) + radiusR) * radiusS radiusR * (sum(W) + radiusS)]) / L;
%     pTest = W * L / ((W + radiusS) * (W + radiusR));

    [A, N] = deal(zeros(1, numEdges));
    for i = 1:numEdges
        A(i) = (sum(W(i:end)) + radiusR) * (radiusS + sum(W(1:i - 1))) / L;
        rS(i) = (radiusS + sum(W(1:i-1))) / L;
        rR(i) = (sum(W(i:end)) + radiusR) / L;
    end
    for i = 1:numEdges -1
        rs = radiusS + sum(W(1:i - 1)); % Segment before w segment
        rr = sum(W(i + 1:end)) + radiusR; % Segement after W segment
        p(i) = W(i) * L / ((W(i) + rs) * (W(i) + rr));
    end

    theta = thetaR - thetaS;
    v = wedgeIndex / 180;

    for i = 1:numEdges
        if theta < 180 - v * 180
            N(i) = -1;
        elseif theta > 180 + v * 180
            N(i) = 1;
        end
    end
    X = A .* cosd((2 .* N .* v .* 180 - theta) / 2) .^ 2;
    
    B = ones(1, numEdges);
    for i = 1:numEdges - 1
        if X(i) <= X(i + 1)
            B(i) = B(i) * p(i);
        else
            B(i + 1) = B(i + 1) * p(i);
        end
    end

    source = [source, zS * ones(numEdges, 1)];
    receiver = [receiver, zR * ones(numEdges, 1)];
end
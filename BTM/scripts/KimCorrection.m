function [scale, A, B] = KimCorrection(data, numEdges, withCorrection)

    radiusS = data.rS;
    W = data.W;
    radiusR = data.rR;
    L = data.L;
    thetaS = data.thetaS;
    thetaR = data.thetaR;
    wedgeIndex = data.wedgeIndex;

    epsilon = 1e-5;
    tS = [thetaS, epsilon * ones(1, numEdges - 1)];
    tR = [wedgeIndex(1:numEdges - 1) - epsilon, thetaR];

    p = zeros(1,numEdges - 1);
    [A, N] = deal(zeros(1, numEdges));
    for i = 1:numEdges
        A(i) = (sum(W(i:end)) + radiusR) * (radiusS + sum(W(1:i - 1))) / L;
    end

    for i = 1:numEdges - 1
        rs = radiusS + sum(W(1:i - 1)); % Segment before w segment
        rr = sum(W(i + 1:end)) + radiusR; % Segement after W segment
        p(i) = W(i) * L / ((W(i) + rs) * (W(i) + rr));
    end
    theta = tR - tS;
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

    if withCorrection
        scale = 1 ./ sqrt(A .* B);
    else
        scale = ones(1, numEdges);
    end
end
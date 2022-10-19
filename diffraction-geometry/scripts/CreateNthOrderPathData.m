function [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, rS, rR, W, height)

    numEdges = length(wedgeIndex);
    W = [mean(W) / 2, W, mean(W) / 2];

    % Create source
    source = [rS * sind(thetaS), rS * cosd(thetaS) - W(1), 0];
    
    % Create edge coordinates
    Q = zeros(numEdges + 1, 3);
    %Q(1,:) = [0, 0.1, 0];
    arg = [0; cumsum(wedgeIndex)];
    for i = 2:numEdges + 2
        vector = W(i - 1) * ([sind(arg(i - 1) + 180 * (i - 3)), cosd(arg(i - 1) + 180 * (i - 3)), 0]);
        Q(i,:) = Q(i - 1,:) + vector;
    end
    
    apex = zeros(numEdges, 3);
    epsilon = 1e-3;
    for i = 1:numEdges
        vector = epsilon * ([sind(arg(i) + wedgeIndex(i) / 2 + 180 * (i - 1)), cosd(arg(i) + wedgeIndex(i) / 2 + 180 * (i - 1)), 0]);
        apex(i,:) = Q(i + 1,:) + vector;
    end

    % Create receiver
    vector = rR * ([sind(arg(numEdges) + 180 * (numEdges - 1) + thetaR), cosd(arg(numEdges) + 180 * (numEdges - 1) + thetaR), 0]);
    receiver = Q(numEdges + 1,:) + vector;
    
    % Check path is valid (no intersections)
    valid = true;
    for i = 1:numEdges - 1
        A = Q(i,1:2);
        B = Q(i + 1,1:2);
        a1 = B(2) - A(2);
        b1 = A(1) - B(1);
        c1 = a1 * A(1) + b1 * A(2);
        pairs = [i, (i + 2:numEdges + 1)];
        for j = 2:length(pairs)
            C = Q(pairs(j),1:2);
            D = Q(pairs(j) + 1,1:2);
            a2 = D(2) - C(2);
            b2 = C(1) - D(1);
            c2 = a2 * C(1) + b2 * C(2);
            determinant = a1*b2 - a2*b1;
            x = (b2*c1 - b1*c2) / determinant;
            y = (a1*c2 - a2*c1) / determinant;
    
            if min(A(1), B(1)) < x && x < max(A(1), B(1)) && min(C(1), D(1)) < x && x < max(C(1), D(1))
                valid = false;
            end
            if min(A(2), B(2)) < y && y < max(A(2), B(2)) && min(C(2), D(2)) < y && y < max(C(2), D(2))
                valid = false;
            end
        end
    end

    % Create CAD file data
    corners = [Q; Q];
    corners(numEdges + 3:end,3) = height;
    
    [source(3), receiver(3), Q(:,3), apex(:,3)] = deal(height / 2);
    
    numEdges = numEdges + 2;
    planeCorners = zeros(numEdges + 2, numEdges);
    % Iterate around the path
    for i = 1:numEdges - 1
        planeCorners(i,1:4) = [i, i + numEdges, i + 1 + numEdges, i + 1];
    end
    % Connect first and last edges
    planeCorners(numEdges, 1:4) = [numEdges, 2 * numEdges, numEdges + 1, 1];
    % Top and bottom plane
    planeCorners(numEdges + 1,:) = 1:numEdges;
    planeCorners(numEdges + 2,:) = [numEdges + 1, fliplr(numEdges + 2:2 * numEdges)];
    
    planeRigid = [ones(1, numEdges - 1), 0, 0, 0];
end
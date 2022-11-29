function [path, direction, heading, data, speed] = CreatePath(landmarks, updateRate, frameRate, scene, wedgeIndex, source, receiverHeading, speed)
    numSegements = size(landmarks,1) - 1;
    walkingSpeed = 1.3;
    if nargin < 8
        speed = walkingSpeed;
    end

    ratio = updateRate / frameRate;
    
    distancePerUpdate = speed / updateRate;
    count = 1;
    path(1,:) = landmarks(1,:);
    vector = landmarks(2,:) - landmarks(1,:);
    direction(1,:) = vector / norm(vector);
    heading(count,:) = receiverHeading;
    data = [];
    for i = 1:numSegements
        vector = landmarks(i + 1,:) - landmarks(i,:);
        distance = norm(landmarks(i + 1,:) - landmarks(i,:));
        stepVector = distancePerUpdate * vector / distance;
        segmentDirection = vector / norm(vector);
        numUpdates = floor(distance / distancePerUpdate);
        for j = 1:numUpdates
            count = count + 1;
            path(count,:) = path(count - 1,:) + stepVector;
            direction(count,:) = segmentDirection;
            heading(count,:) = receiverHeading;
            % data = [data; segmentDirection; receiverHeading];
        end
        landmarks(i + 1,:) = path(count,:);
    end

    for i = 1:count
        idx = ceil((i - 1 +ratio) / ratio);
        saveIdx = 2 * (idx - 1) + 1:2 * idx;
        data(saveIdx,:) = [direction(idx,:); heading(idx,:)];
    end

    start = path(1,:);
    writePath = ['unity-control\', scene, '.csv'];
    writematrix(frameRate, writePath)
    writematrix([wedgeIndex; speed], writePath, 'WriteMode', 'append')
    writematrix([source; receiverHeading; start; data], writePath, 'WriteMode', 'append')
end
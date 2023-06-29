function [azimuth, elevation] = CalculateAzimuthElevation(heading, receiverPosition, sourcePosition)
    
    if size(heading, 1) == 1
        heading = heading .* ones(size(sourcePosition));
    end
    structInput = isstruct(sourcePosition);

    if structInput
        fields = fieldnames(sourcePosition);
        numFields = length(fields);
        for i = 1:numFields
            field = fields{i};
            
            sourceToReceiver = sourcePosition.(field) - receiverPosition;
            sourceAzimuth = CalculateAzimuth(sourceToReceiver(:,1), sourceToReceiver(:,2));
            headingAzimuth = CalculateAzimuth(heading(:,1), heading(:,2));

            idx = sourceAzimuth >= headingAzimuth;
            azimuth.(field) = 360 + sourceAzimuth - headingAzimuth;
            azimuth.(field)(idx) = sourceAzimuth(idx) - headingAzimuth(idx);

            elevation.(field) = CalculateElevation(sourceToReceiver(:,1), sourceToReceiver(:,2), sourceToReceiver(:,3));
        end
    else
        sourceToReceiver = sourcePosition - receiverPosition;
        sourceAzimuth = CalculateAzimuth(sourceToReceiver(:,1), sourceToReceiver(:,2));
        headingAzimuth = CalculateAzimuth(heading(:,1), heading(:,2));

        idx = sourceAzimuth >= headingAzimuth;
        azimuth = 360 + sourceAzimuth - headingAzimuth;
        azimuth(idx) = sourceAzimuth(idx) - headingAzimuth(idx);

        elevation = CalculateElevation(sourceToReceiver(:,1), sourceToReceiver(:,2), sourceToReceiver(:,3));
    end
end
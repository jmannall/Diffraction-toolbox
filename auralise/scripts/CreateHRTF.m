function ir = CreateHRTF(azimuth, elevation)

    load 'ReferenceHRTF.mat' hrtfData sourcePosition

    hrtfData = permute(double(hrtfData),[2,3,1]);
    
    sourcePosition = sourcePosition(:,[1,2]);

    structInput = isstruct(azimuth);

    if structInput
        fields = fieldnames(azimuth);
        numFields = length(fields);
        for i = 1:numFields
            idx = fields{i};
            targetPosition = [azimuth.(idx) elevation.(idx)];

            interpolatedIR  = interpolateHRTF(hrtfData,sourcePosition,targetPosition);

            ir.(idx).L = squeeze(interpolatedIR(:,1,:))';
            ir.(idx).R = squeeze(interpolatedIR(:,2,:))';
        end
    else
        targetPosition = [azimuth elevation];
        
        interpolatedIR  = interpolateHRTF(hrtfData,sourcePosition,targetPosition);
        
        ir.L = squeeze(interpolatedIR(:,1,:))';
        ir.R = squeeze(interpolatedIR(:,2,:))';
    end
end
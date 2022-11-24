%% Create CAD file from corners and planeCorners

function cadFilePath = CreateCADFile(inFilePath, index, corners, planeCorners, planeRigid)

    CheckFileDir('geometry')
    
    cadFilePath = [inFilePath, filesep, 'geometry', filesep, num2str(index), '_geo.cad'];

    % Open CAD file
    cadfile = fopen(cadFilePath, 'w');

    % Write CAD file
    numCorners = length(corners);
    heading = '%CORNERS';
    fprintf(cadfile, '%s\n', heading);
    count = 1;
    for i = 1:numCorners
        fprintf(cadfile, '%d %2.4f %2.4f %2.4f\n', count, corners(i,:));
        count = count + 1;
    end
    
    [numPlanes, cornersPerPlane] = size(planeCorners);
    heading = '%PLANES';
    fprintf(cadfile, '\n%s\n', heading);
    count = 1;
    for i = 1:numPlanes
        if(planeRigid(i) == 0)
            line = ' / /TOTABS';
        else
            line = ' / /RIGID';
        end
        inputs = blanks(3 * cornersPerPlane);
        for j = 1:cornersPerPlane
            inputs(3 * (j - 1) + 1:3 * j) = ' %d';
        end
        fprintf(cadfile, [' %d %s\n', inputs, '\n \n'], count, line, planeCorners(i,:));
        count = count + 1;
    end

    % Close CAD file
    fclose(cadfile);
end
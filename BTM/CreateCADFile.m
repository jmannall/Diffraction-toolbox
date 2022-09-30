%% Create CAD file from corners and planeCorners

function cadFilePath = CreateCADFile(inFilePath, index, corners, planeCorners, planeRigid)

    if ~exist('geometry', 'dir')
       mkdir geometry
    end
    
    cadFilePath = [inFilePath, '\geometry\', num2str(index), '_geo.cad'];

    % Open CAD file
    cadfile = fopen(cadFilePath, 'w');

    % Write CAD file
    heading = '%CORNERS';
    fprintf(cadfile, '%s\n', heading);
    count = 1;
    for i = 1:length(corners)
        fprintf(cadfile, '%d %2.4f %2.4f %2.4f\n', count, corners(i,:));
        count = count + 1;
    end
    
    heading = '%PLANES';
    fprintf(cadfile, '\n%s\n', heading);
    count = 1;
    for i = 1:size(planeCorners,1)
        if(planeRigid(i) == 0)
            line = ' / /TOTABS';
        else
            line = ' / /RIGID';
        end
        fprintf(cadfile, ' %d %s\n %d %d %d %d\n \n', count, line, planeCorners(i,:));
        count = count + 1;
    end

    % Close CAD file
    fclose(cadfile);
end
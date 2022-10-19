% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create ir, tf data for a single wall in free field

% All coordinates go x, y, z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ir, tfmag, tvec, fvec, tfcomplex] = SingleWall(wallHeight,wallThickness,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters,createPlot)
    
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({wallHeight,wallThickness,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters});
    [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(mFile, index);

    if resultExists
        load(loadPath, "ir", "tfmag", "tvec", "fvec", "tfcomplex");
        disp('IR load from save');
    else
        
        % Check for invalid data
        wallIndex = 360 - asind((wallThickness / 2) / radiusR);
        if thetaR >= wallIndex
            error('Receiver angle exceeds the exterior angle of the wall');
        end
        if thetaS <= (360 - wallIndex) / 2
            error('Sources lies within the bounds of the wall');
        end

        % Create geometry data
        wallSize = 10 * max(radiusS, radiusR);
        corners = [0 0 0
            0 -wallThickness 0
            wallSize 0 0
            wallSize -wallThickness 0
            0 0 wallHeight
            0 -wallThickness wallHeight
            wallSize 0 wallHeight
            wallSize -wallThickness wallHeight];

        planeCorners = [1 3 4 2
            1 2 6 5
            1 5 7 3
            5 6 8 7
            3 7 8 4
            2 4 8 6];

        planeRigid = [0 1 1 0 0 1];

%         source = [radiusS * cosd(thetaS) + wallThickness / 2, radiusS * sind(thetaS) - wallThickness / 2, zS];
%         receiver = [radiusR * cosd(thetaR) + wallThickness / 2, radiusR * sind(thetaR) - wallThickness / 2, zR];
        source = [radiusS * cosd(thetaS) + wallThickness / 2, radiusS * sind(thetaS) - wallThickness / 2, zS];
        receiver = [radiusR * cosd(thetaR) + wallThickness / 2, radiusR * sind(thetaR) - wallThickness / 2, zR];

        % Plot geometry
        if createPlot
            PlotGeometry(corners, planeCorners, source, receiver)
        end
        
        % Create CAD file   
        cadFilePath = CreateCADFile(inFilePath, index, corners, planeCorners, planeRigid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Find BTM response
        geofiledata = struct('geoinputfile',cadFilePath);
        Sindata = struct('coordinates',source);
        Rindata = struct('coordinates',receiver);
        controlparameters.savealldifforders = 1;
        filehandlingparameters = struct('outputdirectory',[inFilePath,filesep,'results']);
        filehandlingparameters.filestem = [fileName, '_', index];
        filehandlingparameters.savelogfile = 0;
        filehandlingparameters.showtext = 1;
        filehandlingparameters.suppressresultrecycling = 1;
    
        EDmain_convex_time(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %     controlparameters.difforder = 15;
    %     controlparameters.ngauss = 24;
    %     fvec_tf = fvec(1:100);
    %     controlparameters.frequencies = fvec_tf;
    %     EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Process and save the results
        [ir, tfmag, tfcomplex, tvec, fvec] = ProcessBTMResults(inFilePath, filehandlingparameters, controlparameters, cadFilePath, savePath);
    end
end
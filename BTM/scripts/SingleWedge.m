% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create ir, tf data for a single wedge in free field

% All coordinates go x, y, z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters,createPlot)

    % Create file info
    mFile = mfilename('fullpath');
    controlparameters.difforder = 1;
    index = DataHash({wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters});
    [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(mFile, index);

    if resultExists
        load(loadPath, "ir", "tfmag", "tvec", "fvec", "tfcomplex");
        disp('IR load from save');
    else
        
        % Check for invalid data
        if wedgeIndex > 360
            error('Wedge index exceeds 360 degrees');
        elseif wedgeIndex <= 180
            error('Wedge is concave');
        end
        if thetaR >= wedgeIndex
            error('Receiver angle exceeds the wedge index');
        end
        if thetaS == 0
            error('Source angle lies on a plane');
        end
                
        % Create geometry data
        wedgeSize = 100 * max(max(radiusS), max(radiusR));
        x = wedgeSize * cosd(wedgeIndex);
        y = wedgeSize * sind(wedgeIndex);
        corners = [0 0 0
            wedgeSize 0 0
            x y 0
            0 0 wedgeLength
            wedgeSize 0 wedgeLength
            x y wedgeLength
            wedgeSize + 1 y 0
            wedgeSize + 1 y wedgeLength];
        
        planeCorners = [1 3 6 4
            3 7 8 6
            2 1 4 5
            4 6 8 5
            1 2 7 3
            2 5 8 7];

        planeRigid = [1 0 1 0 0 0];
        
        source = [radiusS * cosd(thetaS), radiusS * sind(thetaS), zS];
        receiver = [radiusR .* cosd(thetaR), radiusR .* sind(thetaR), zR];

        % Plot geometry 
        if createPlot
            if isfield(controlparameters, 'Rstart')
                receiverPlot = receiver / controlparameters.Rstart;
            else
                receiverPlot = receiver;
            end
            PlotGeometry(corners, planeCorners, source, receiverPlot)
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
        filehandlingparameters.showtext = 0;
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
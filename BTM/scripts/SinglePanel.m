% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create ir, tf data for a single panel in free field

% All coordinates go x, y, z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ir, tfmag, tvec, fvec, tfcomplex] = SinglePanel(width, depth, height, controlparameters, createPlot)
    
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({width,depth,height,controlparameters});
    [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(mFile, index);

    if resultExists
        load(loadPath, "ir", "tfmag", "tvec", "fvec", "tfcomplex");
        disp('IR load from save');
    else

        % Create geometry data
        x = width / 2;
        corners = [-x 0 0       %1
            x 0 0               %2
            -x -depth 0         %3  
            x -depth 0          %4
            -x 0 height         %5
            x 0 height          %6
            -x -depth height    %7
            x -depth height];   %8
        
        planeCorners = [1 2 4 3
            1 5 6 2
            1 3 7 5
            5 7 8 6
            3 4 8 7
            2 6 8 4];
    
        planeRigid = [1 1 1 1 1 1];
        
        rS = 2.5;
        thetaS = -45;
        source = [rS * sind(thetaS), rS * cosd(thetaS), height / 2 * ones(size(thetaS))];

        rR = 0.8;
        thetaR = (0.1:1:180.1)';
        receiver = [rR * sind(thetaR), rR * cosd(thetaR), height / 2 * ones(size(thetaR))];
    
        % Plot geometry
        if createPlot
            PlotGeometry(corners, planeCorners, source, receiver)
        end    
        % Create CAD file
        cadFilePath = CreateCADFile(inFilePath, index, corners, planeCorners, planeRigid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Find BTM response
        geofiledata = struct('geoinputfile',cadFilePath);
        Sindata = struct('coordinates',source,'noDirect',controlparameters.noDirect);
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
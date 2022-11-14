% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create ir, tf data for a single wedge in free field

% All coordinates go x, y, z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ir, tfmag, tvec, fvec, tfcomplex] = DefaultBTM(controlparameters)

    wedgeLength = 10;
    wedgeIndex = 270;
    thetaS = 10;
    thetaR = 250;
    radiusS = 1;
    radiusR = 1;
    zS = 5;
    zR = 5;

    if ~exist('results\SingleWedge', 'dir')
           mkdir results\SingleWedge
    end

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
        wedgeSize = 10 * max(radiusS, radiusR);
        corners = [0 0 0
            wedgeSize 0 0
            wedgeSize * cosd(wedgeIndex) wedgeSize * sind(wedgeIndex) 0
            0 0 wedgeLength
            wedgeSize 0 wedgeLength
            wedgeSize * cosd(wedgeIndex) wedgeSize * sind(wedgeIndex) wedgeLength
            wedgeSize -0.0001 0
            wedgeSize -0.0001 wedgeLength];
        
        planeCorners = [1 3 6 4
            3 7 8 6
            2 1 4 5
            4 6 8 5
            1 2 7 3
            2 5 8 7];

        planeRigid = [1 0 1 0 0 0];
        
        source = [radiusS * cosd(thetaS), radiusS * sind(thetaS), zS];
        receiver = [radiusR * cosd(thetaR), radiusR * sind(thetaR), zR];
    
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
        ir = [ir.diff1];
        tfmag = [tfmag.diff1];
        tfcomplex = [tfcomplex.diff1];
    end
end
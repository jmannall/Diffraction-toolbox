function [ir, tfmag, tvec, fvec, tfcomplex] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot)

    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({source, receiver, corners, planeCorners, planeRigid});
    [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(mFile, index);
    
    if resultExists
        load([cd filesep loadPath], "ir", "tfmag", "tvec", "fvec", "tfcomplex");
        disp('IR load from save');
    else

        % Plot geometry
        if createPlot
            PlotGeometry(corners, planeCorners, source, receiver)
        end

        % Create CAD file
        cadFilePath = CreateCADFile(inFilePath, index, corners, planeCorners, planeRigid);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Find BTM response
        geofiledata = struct('geoinputfile',cadFilePath);
        %geofiledata = struct('corners',corners, 'planecorners', planeCorners);
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

        % Process and save the results
        [ir, tfmag, tfcomplex, tvec, fvec] = ProcessBTMResults(inFilePath, filehandlingparameters, controlparameters, cadFilePath, savePath);
    end
end
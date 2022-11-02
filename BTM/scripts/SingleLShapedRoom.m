function [tfmag, fvec, tfcomplex] = SingleLShapedRoom(x, y, height, receiver, controlparameters, outRoom, createPlot)

    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({x, y, height, receiver, controlparameters});
    [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(mFile, index);
    
    if outRoom
        n = controlparameters.nfft / 2;
        template.complete = zeros(n, 1);
        template.direct = zeros(n, 1);
        template.geom = zeros(n, 1);
        difforder = controlparameters.difforder;
        for i = 1:difforder
            idx = ['diff', num2str(i)];
            template.(idx) = zeros(n, 1);
        end
        [tfmag, tfcomplex] = deal(repmat(template, 1, 1));
        fvec = controlparameters.fs/controlparameters.nfft*[0:controlparameters.nfft/2-1];
    else
        if resultExists
            load(loadPath, "ir", "tfmag", "tvec", "fvec", "tfcomplex");
            disp('IR load from save');
        else
                [corners, planeCorners, source, receiver] = CreateLShapedRoomGeometry(x, y, height, receiver);
        
                % Plot geometry
                if createPlot
                    PlotGeometry(corners, planeCorners, source, receiver)
                end
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                % Find BTM response
                geofiledata = struct('corners',corners, 'planecorners', planeCorners);
                Sindata = struct('coordinates',source);
                Rindata = struct('coordinates',receiver);
                controlparameters.savealldifforders = 1;
                filehandlingparameters = struct('outputdirectory',[inFilePath,filesep,'results']);
                filehandlingparameters.filestem = [fileName, '_', index];
                filehandlingparameters.savelogfile = 0;
                filehandlingparameters.showtext = 1;
                filehandlingparameters.suppressresultrecycling = 1;
            
                EDmain_nonconvex_time(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                % Process and save the results
                cadFilePath = 'none';
                [ir, tfmag, tfcomplex, tvec, fvec] = ProcessBTMResults(inFilePath, filehandlingparameters, controlparameters, cadFilePath, savePath);
        end
    end
end
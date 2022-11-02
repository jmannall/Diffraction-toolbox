function [tfmag, fvec, tfcomplex] = SingleDiffractionGrating(numPillars, gratingWidth, pillarWidth, pillarHeight, receiver, controlparameters, createPlot)

    receiver(2) = abs(receiver(2));
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({numPillars, gratingWidth, pillarWidth, pillarHeight, receiver, controlparameters});
    [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(mFile, index);
    
    if resultExists
        load(loadPath, "tfmag", "fvec", "tfcomplex");
        disp('IR load from save');
    else

        [corners, planeCorners, source, receiver] = CreateDiffractionGratingGeometry(numPillars, gratingWidth, pillarWidth, pillarHeight, receiver);

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
        controlparameters.ngauss = 24;
        nfft = controlparameters.nfft;
        fs = controlparameters.fs;
        fvec = fs/nfft*[0:nfft/2-1];
        fvec = [84 369 1625];
        controlparameters.frequencies = fvec;
        filehandlingparameters = struct('outputdirectory',[inFilePath,filesep,'results']);
        filehandlingparameters.filestem = [fileName, '_', index];
        filehandlingparameters.savelogfile = 0;
        filehandlingparameters.showtext = 1;
        filehandlingparameters.suppressresultrecycling = 1;
        
        EDmain_nonconvexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Process and save the results
        % [ir, tfmag, tfcomplex, tvec, fvec] = ProcessBTMResults(inFilePath, filehandlingparameters, controlparameters, cadFilePath, savePath);
    
        load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'], 'tfinteqdiff');
        load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'], 'tfdirect', 'tfgeom', 'tfdiff');

        % Create template
        template.complete = [];
        template.direct = [];
        template.geom = [];
        template.diff1 = [];
        template.diffHod = [];
        [tfmag, tfcomplex] = deal(repmat(template, 1, 1));
        tfcomplex.complete = tfdirect + tfgeom + tfdiff + tfinteqdiff;
        tfcomplex.direct = tfdirect;
        tfcomplex.geom = tfgeom;
        tfcomplex.diff1 = tfdiff;
        tfcomplex.diffHod = tfinteqdiff;

        tfmag.complete = mag2db(abs(tfcomplex.complete));
        tfmag.direct = mag2db(abs(tfcomplex.direct));
        tfmag.geom = mag2db(abs(tfcomplex.geom));
        tfmag.diff1 = mag2db(abs(tfcomplex.diff1));
        tfmag.diffHod = mag2db(abs(tfcomplex.diffHod));
                
        save(savePath, "tfmag", "fvec", "tfcomplex");
        disp('Result saved')
        
        % Clear up and delete files
        %delete(cadFilePath);
        path = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem];
        delete([path, '_eddata.mat']);
        delete([path, '_ed2data.mat']);
        delete([path, '_submatrixdata.mat']);
        delete([path, '_paths.mat']);
        delete([path, '_Rdata.mat']);
        delete([path, '_Sdata.mat']);
        delete([path, '_tfinteq.mat']);
        delete([path, '_tf.mat']);
    end
end
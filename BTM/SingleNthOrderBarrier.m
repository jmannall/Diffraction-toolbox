% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create ir, tf data for a single Nth order barrier in free field

% All coordinates go x, y, z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ir, tfmag, tvec, fvec, tfcomplex] = SingleNthOrderBarrier(barrierRadius,barrierHeight,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters,createPlot)
    
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({barrierHeight,barrierRadius,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters});
    [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(mFile, index);

    if resultExists
        load(loadPath, "ir", "tfmag", "tvec", "fvec", "tfcomplex");
        disp('IR load from save');
    else
        
        % Create geometry data
        numDiffEdges = controlparameters.difforder;
        
        [corners, planeCorners, planeRigid, source, receiver] = CreateNthOrderBarrierGeometry(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, numDiffEdges);
       
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
        
%         controlparameters.ngauss = 24;
%         nfft = controlparameters.nfft;
%         fs = controlparameters.fs;
%         fvec = fs/nfft*[0:nfft/2-1];
%         fvec_tf = fvec(1:500);
%         controlparameters.frequencies = fvec_tf;
%         EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'], 'tfinteqdiff');
%         
%         tfcomplex = tfinteqdiff;
%         tfmag = mag2db(abs(tfcomplex));
%         ir = 0;
%         tvec = 0;
% 
%         save(savePath, "ir", "tfmag", "tvec", "fvec_tf", "tfcomplex");
%         disp('Result saved')
%     
%         % Clear up and delete files
%         delete(cadFilePath);
%         path = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem];
%         delete([path, '_eddata.mat']);
%         delete([path, '_ed2data.mat']);
%         delete([path, '_submatrixdata.mat']);
%         delete([path, '_ir.mat']);
%         delete([path, '_paths.mat']);
%         delete([path, '_Rdata.mat']);
%         delete([path, '_Sdata.mat']);
%         delete([path, '_tfinteq.mat']);
%         delete([path, '_tf.mat']);
        % Process and save the results
        [ir, tfmag, tfcomplex, tvec, fvec] = ProcessBTMResults(inFilePath, filehandlingparameters, controlparameters, cadFilePath, savePath);
    end
end
% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create ir, tf data for a single wedge in free field

% All coordinates go x, y, z

function [ir, tf, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('results\SingleWedge', 'dir')
           mkdir results\SingleWedge
    end

    % Create file info
    mFile = mfilename('fullpath');
    [inFilePath,fileName] = fileparts(mFile);
    index = DataHash([wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs]);
    fileStem = [fileName, '_', num2str(index)];
    savepath = ['results\SingleWedge\', fileStem];
    loadpath = ['results\SingleWedge\', fileStem, '.mat'];

    test = exist(loadpath, "file");

    if test == 2
        load(loadpath, "ir", "tf", "tvec", "fvec", "tfcomplex");
        disp('IR load from save');
        save(savepath, "ir", "tf", "tvec", "fvec", "tfcomplex");
    else
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Input data
        wedgeSize = 10 * max(radiusS, radiusR);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Check for invalid data
        if wedgeIndex > 360
            disp('Wedge index exceeds 360 degrees');
            return
        elseif wedgeIndex <= 180
            disp ('Wedge is concave');
            return
        end
        
        if thetaR >= wedgeIndex
            disp('Receiver angle exceeds the wedge index');
            return
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Create geometry data
        corners = [0 0 0
            wedgeSize 0 0
            wedgeSize * cosd(wedgeIndex) wedgeSize * sind(wedgeIndex) 0
            0 0 wedgeLength
            wedgeSize 0 wedgeLength
            wedgeSize * cosd(wedgeIndex) wedgeSize * sind(wedgeIndex) wedgeLength
            wedgeSize -0.0001 0
            wedgeSize -0.0001 wedgeLength];
        
        planecorners = [1 3 6 4
            3 7 8 6
            2 1 4 5
            4 6 8 5
            1 2 7 3
            2 5 8 7];

        planerigid = [1 0 1 0 0 0];
        
        source = [radiusS * cosd(thetaS), radiusS * sind(thetaS), zS];
        receiver = [radiusR * cosd(thetaR), radiusR * sind(thetaR), zR];

        plane = corners(planecorners',:);
        
%         figure(1)
%         plot3(source(1), source(2), source(3), 'o');
%         hold on
%         plot3(receiver(1), receiver(2), receiver(3), 'o');
%         legend('Source', 'Receiver')
%         for i = 1:4:length(plane)
%             planePlot = [plane(i:i + 3,:); plane(i,:)];
%             plot3(planePlot(:,1), planePlot(:,2), planePlot(:,3))
%         end
%         xlim([-1 1])
%         ylim([-1 1])
%         grid on
%         hold off
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create CAD file
        if ~exist('geometry', 'dir')
           mkdir geometry
        end
        
        csvFilePath = [inFilePath, filesep, 'geometry', filesep];
        
        map = num2str(index);
        cadFilePath = [csvFilePath, map, '_geo.cad'];
        
        cadfile = fopen(cadFilePath, 'w');

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
        for i = 1:size(planecorners,1)
            if(planerigid(i) == 0)
                line = ' / /TOTABS';
            else
                line = ' / /RIGID';
            end
            fprintf(cadfile, ' %d %s\n %d %d %d %d\n \n', count, line, planecorners(i,:));
            count = count + 1;
        end
        fclose(cadfile);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Find BTM response
        geofiledata = struct('geoinputfile',cadFilePath);
        Sindata = struct('coordinates',source);
        Rindata = struct('coordinates',receiver);
        controlparameters = struct('fs',fs);
        controlparameters.difforder = 1;
        controlparameters.savealldifforders = 1;
        filehandlingparameters = struct('outputdirectory',[inFilePath,filesep,'results']);
        filehandlingparameters.filestem = [fileName, index];
        filehandlingparameters.savelogfile = 0;
        filehandlingparameters.showtext = 1;
        filehandlingparameters.suppressresultrecycling = 1;

        EDmain_convex_time(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
        
        nfft = 8192;
        fvec = controlparameters.fs/nfft*[0:nfft/2-1];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %     controlparameters.difforder = 1;
    %     controlparameters.ngauss = 24;
    %     fvec_tf = fvec(1:100);
    %     controlparameters.frequencies = fvec_tf;
    %     EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Load and save the results

        %eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_ir.mat'''])
        load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_ir.mat'], 'irdirect', 'irgeom', 'irdiff');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Present the results
        template.complete = [];
        template.direct = [];
        template.geom = [];
        template.diff = [];

        [ir, tf, tfcomplex] = deal(repmat(template, 1, 1));

        ir.complete = irdirect + irgeom + irdiff;
        ir.direct = irdirect;
        ir.geom = irgeom;
        ir.diff = irdiff;

        [tf.complete, tfcomplex.complete] = irToTf(ir.complete, nfft);
        [tf.direct, tfcomplex.direct] = irToTf(irdirect, nfft);
        [tf.geom, tfcomplex.geom] = irToTf(irgeom, nfft);
        [tf.diff, tfcomplex.diff] = irToTf(irdiff, nfft);

        ndiff = length(ir.complete);
        tvec = 1/controlparameters.fs*[0:ndiff-1];

        save(savepath, "ir", "tf", "tvec", "fvec", "tfcomplex");
        delete(cadFilePath);
        path = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem];
        delete([path, '_eddata.mat']);
        delete([path, '_ir.mat']);
        delete([path, '_paths.mat']);
        delete([path, '_Rdata.mat']);
        delete([path, '_Sdata.mat']);
    end
end
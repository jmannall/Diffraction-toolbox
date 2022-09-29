% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create ir, tf data for a single wall in free field

% All coordinates go x, y, z

function [ir, tf, tvec, fvec, tfcomplex] = SingleWall(wallHeight,wallThickness,thetaS,thetaR,radiusS,radiusR,zS,zR,fs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create file info
    mFile = mfilename('fullpath');
    [inFilePath,fileName] = fileparts(mFile);
    index = DataHash([wallHeight,wallThickness,thetaS,thetaR,radiusS,radiusR,zS,zR,fs]);
    fileStem = [fileName, '_', num2str(index)];
    savepath = ['results\SingleWall\', fileStem];
    loadpath = ['results\SingleWall\', fileStem, '.mat'];

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
        
        % Create geometry data
        corners = [0 0 0
            0 -wallThickness 0
            wedgeSize 0 0
            wedgeSize -wallThickness 0
            0 0 wallHeight
            0 -wallThickness wallHeight
            wedgeSize 0 wallHeight
            wedgeSize -wallThickness wallHeight];

        planecorners = [1 3 4 2
            1 2 6 5
            1 5 7 3
            5 6 8 7
            3 7 8 4
            2 4 8 6];

        planerigid = [0 1 1 0 0 1];
        
        radiusS = radiusS + wallThickness / 2;
        radiusR = radiusR + wallThickness / 2;
        source = [radiusS * cosd(thetaS) + wallThickness / 2, radiusS * sind(thetaS) - wallThickness / 2, zS];
        receiver = [radiusR * cosd(thetaR) + wallThickness / 2, radiusR * sind(thetaR) - wallThickness / 2, zR];
        plane = corners(planecorners',:);
        
%         figure(1)
%         plot3(source(1), source(2), source(3), 'o');
%         hold on
%         plot3(receiver(1), receiver(2), receiver(3), 'o');
%         legend('Source', 'Receiver','source1', 'receiver1')
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
        controlparameters.difforder = 3;
        controlparameters.savealldifforders = 1;
        filehandlingparameters = struct('outputdirectory',[inFilePath,filesep,'results']);
        filehandlingparameters.filestem = [fileName, index];
        filehandlingparameters.savelogfile = 0;
        filehandlingparameters.showtext = 1;
        filehandlingparameters.suppressresultrecycling = 0;

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
        load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_irhod.mat'], 'irhod');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        irdiffTwo = irhod{2};
        irdiffThree = irhod{3};

        nir = length(irdirect);
        ndiff = length(irdiffTwo);

        if ndiff > nir
            padding = zeros(ndiff - nir, 1);
            irdirect = [irdirect; padding];
            irgeom = [irgeom; padding];
            irdiff = [irdiff; padding];
        else
            irdiffTwo = zeros(size(irdirect));
            irdiffThree = zeros(size(irdirect));
        end

        % Present the results
        template.complete = [];
        template.direct = [];
        template.geom = [];
        template.diff = [];
        template.diffTwo = [];
        template.diffTree = [];

        [ir, tf, tfcomplex] = deal(repmat(template, 1, 1));

        [ir.complete, tf.complete, tfcomplex.complete] = ComputeWedgeResponses(irdirect + irgeom + irdiff + irdiffTwo + irdiffThree, nfft);
%         [ir.complete, tf.complete, tfcomplex.complete] = ComputeWedgeResponses(irdirect + irgeom + irdiff + irdiffTwo, nfft);
        [ir.direct, tf.direct, tfcomplex.direct] = ComputeWedgeResponses(irdirect, nfft);
        [ir.geom, tf.geom, tfcomplex.geom] = ComputeWedgeResponses(irgeom, nfft);
        [ir.diff, tf.diff, tfcomplex.diff] = ComputeWedgeResponses(irdiff, nfft);
        [ir.diffTwo, tf.diffTwo, tfcomplex.diffTwo] = ComputeWedgeResponses(irdiffTwo, nfft);
        [ir.diffThree, tf.diffThree, tfcomplex.diffThree] = ComputeWedgeResponses(irdiffThree, nfft);

        tvec = 1/controlparameters.fs*[0:max(nir,ndiff)-1];

        save(savepath, "ir", "tf", "tvec", "fvec", "tfcomplex");
        delete(cadFilePath);
        path = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem];
        delete([path, '_eddata.mat']);
        delete([path, '_ir.mat']);
        delete([path, '_irhod.mat']);
        delete([path, '_paths.mat']);
        delete([path, '_Rdata.mat']);
        delete([path, '_Sdata.mat']);
    end
end
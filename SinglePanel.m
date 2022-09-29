% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create geometry data for a single wedge in free field

% All coordinates go x, y, z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ir, tf, tvec, fvec, tfcomplex] = SinglePanel(width,depth,height,fs)

% Create file info
mFile = mfilename('fullpath');
[inFilePath,fileName] = fileparts(mFile);
index = DataHash([width,depth,height,fs]);
fileStem = [fileName, '_', num2str(index)];
savepath = ['results\SinglePanel\', fileStem];
loadpath = ['results\SinglePanel\', fileStem, '.mat'];

test = exist(loadpath, "file");

if test == 2
    load(loadpath, "ir", "tf", "tvec", "fvec", "tfcomplex");
    disp('IR load from save');
    save(savepath, "ir", "tf", "tvec", "fvec", "tfcomplex");
else
    
    % Create geometry data
    corners = [-width /2 0 0    %1
        width / 2 0 0           %2
        -width /2 -depth 0      %3  
        width / 2 -depth 0      %4
        -width /2 0 height      %5
        width / 2 0 height      %6
        -width /2 -depth height %7
        width / 2 -depth height];%8
    
    planecorners = [1 2 4 3
        1 5 6 2
        1 3 7 5
        5 7 8 6
        3 4 8 7
        2 6 8 4];

    planerigid = [1 1 1 1 1 1];
    
    source = [0, 1, height / 2];
    receiver = [0, -1, height / 2];

    plane = corners(planecorners',:);
    idx1 = 1:4:24;
    idx2 = 2:4:24;
    idx3 = 3:4:24;
    norms = cross(plane(idx1,:) - plane(idx3,:), plane(idx1,:) - plane(idx2,:));

    norms = 0.1 * normalize(norms);
    figure(1)
    plot3(source(1), source(2), source(3), 'o');
    hold on
    plot3(receiver(1), receiver(2), receiver(3), 'o');
    legend('Source', 'Receiver')
    for j = 1:length(plane) / 4
        i = (j - 1) * 4 + 1;
        %planePlot = [plane(i:i + 3,:); plane(i,:)];
        planePlot = [plane(i:i + 3,:); plane(i,:)];
        planeCentre = [sum(planePlot(1:4,1)) / 4, sum(planePlot(1:4,2)) / 4, sum(planePlot(1:4,3)) / 4];
        normPlot = [planeCentre; planeCentre - norms(j,:)];
        plot3(planePlot(:,1), planePlot(:,2), planePlot(:,3))
        plot3(normPlot(:,1), normPlot(:,2), normPlot(:,3), 'k')
    end
    xlim([-0.5 0.5])
    ylim([-1 1])
    zlim([-0.25 0.75])
    grid on
    hold off

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
    load([inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_ir.mat'],'irdiff', 'irdirect');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Present the results
    
    ir = irdiff;
    ndiff = length(irdiff);
    tvec = 1/controlparameters.fs*[0:ndiff-1];
    
    % eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'''])
    % eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])
    
    F = fft(irdiff,nfft);
    tfcomplex = F(1:nfft/2,:);
    tf = 20*log10(abs(F(1:nfft/2,:)));

    save(savepath, "ir", "tf", "tvec", "fvec", "tfcomplex");
    delete(cadFilePath);
    path = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem];
    delete([path, '_eddata.mat']);
    delete([path, '_ir.mat']);
    delete([path, '_paths.mat']);
    delete([path, '_Rdata.mat']);
    delete([path, '_Sdata.mat']);
end

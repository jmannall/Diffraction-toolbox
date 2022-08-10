% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create geometry data for a single wedge in free field

% All coordinates go x, y, z

function [ir, tf, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create file info
    mFile = mfilename('fullpath');
    [inFilePath,fileStem] = fileparts(mFile);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Input data
    thetaW = 360 - wedgeIndex;
    wedgeSize = max(radiusS, radiusR);
    
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
        wedgeSize * cosd(wedgeIndex) wedgeSize * sind(wedgeIndex) wedgeLength];
    
    planecorners = [1 3 6 4
        3 2 5 6
        2 1 4 5
        4 6 5 0
        1 2 3 0];
    
    planerigid = [1 0 1 0 0];
    
    source = [radiusS * cosd(thetaS), radiusS * sind(thetaS), zS];
    receiver = [radiusR * cosd(thetaR), radiusR * sind(thetaR), zR];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create CAD file
    if ~exist('geometry', 'dir')
       mkdir geometry
    end
    
    csvFilePath = [inFilePath, filesep, 'geometry', filesep];
    
    map = 'test';
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
    filehandlingparameters.filestem = fileStem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 1;
    
    EDmain_convex_time(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
    
    nfft = 4096;
    fvec = controlparameters.fs/nfft*[0:nfft/2-1];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    controlparameters.difforder = 1;
    controlparameters.ngauss = 24;
    fvec_tf = fvec(1:100);
    controlparameters.frequencies = fvec_tf;
    EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load and present the results
    
    eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_ir.mat'''])
    
    ir = irdiff;
    ndiff = length(irdiff);
    tvec = 1/controlparameters.fs*[0:ndiff-1];
    
    eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'''])
    eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])
    
    F = fft(irdiff,nfft);
    tfcomplex = F(1:nfft/2,:);
    tf = 20*log10(abs(F(1:nfft/2,:)));
end
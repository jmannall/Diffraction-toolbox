

% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create geometry data for a single wedge in free field

% All coordinates go x, y, z

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create file info
mFile = mfilename('fullpath');
[inFilePath,fileStem] = fileparts(mFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input data
wedgeLength = 2;
wedgeIndex = 270;
thetaW = 360 - wedgeIndex;
thetaS = 30;
thetaR = 250;
radiusS = 1;
radiusR = 0.8;
zS = 2.2;
zR = 3;
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
controlparameters = struct('fs',48000);
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
    
nreceivers = size(Rindata.coordinates,1);
nsources = size(Sindata.coordinates,1);

eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_ir.mat'''])
%eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_irhod.mat'''])

%ncells = size(irhod,2);
%ntotlength = size(irhod{2},1);
%allirtots = zeros(ntotlength,ncells+1);
%for ii = 2:ncells    
%   allirtots(:,ii+1) = irhod{ii}; 
%end

ndir = length(irdirect);
allirtots(1:ndir,1) = irdirect + irgeom;
allirtots(1:ndir,2) = irdiff;

cumsumirtots = cumsum(allirtots.').';

eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'''])
eval(['load ''',inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])

tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;

%tvec = 1/controlparameters.fs*[0:ntotlength-1];
tvec = 1/controlparameters.fs*[0:ndir-1];

% Plot the model
figure(1)
eddatafile = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
EDplotmodel(eddatafile,3)

% Plot the irs
figure(2)
h = plot(tvec*1e3,cumsumirtots,'-');
% set(h(1),'MarkerSize',2);
set(h(2),'LineWidth',2);
g = get(h(1),'Parent');
set(g,'FontSize',14)
grid
g = xlabel('Time   [ms]');
set(g,'FontSize',14)
g = ylabel('Impulse response   [-]');
set(g,'FontSize',14)
axis([6 12 -0.02 0.003])
g = legend('GA','1st-order diffr.');
set(g,'FontSize',14,'Location','SouthEast')


F = fft(cumsumirtots,nfft);

% Plot the frequency response
figure(3)
h = semilogx(fvec,20*log10(abs(F(1:nfft/2,:))));
for ii = 1:2
   set(h(ii),'LineWidth',2) 
end
g = get(h(1),'Parent');
set(g,'FontSize',14)
grid
g = xlabel('Frequency   [Hz]');
set(g,'FontSize',14)
g = ylabel('Frequency response magnitude   [dB]');
set(g,'FontSize',14)
g = legend('GA only','Incl. diffr.1');
set(g,'FontSize',14,'Location','SouthEast')
axis([20 20000 -8 4])

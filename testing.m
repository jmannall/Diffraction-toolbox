
% Ask Enzo about adding automatic testing features. e.g Run a document and
% return a 1 if everything running as expected.

% Look into adding a datahash for each calculation. Could include in the
% single wedge function so if gets the same wedge dimensions
% (non-consecutively) can still recycle results.

% Split into functions e.g make CAD file.

% Create input struct.

% Requires EDtoolbox by upsvensson, DataHash.m and lgwt.m to run

% Create geometry data for a single wedge in free field

% All coordinates go x, y, z

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;

step = 10;
minw = 180;
maxw = 360;
shadowZone = true;

wedgeIndex = 220;
thetaS = 10;
thetaR = 30;

thetaS = 190;
thetaR = 210;

wedge = 190:10:360;
bendingAngle = 190:10:360;
minAngle = 0:10:180;

result = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, wedge, bendingAngle, minAngle);


test = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, fs);

test2 = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, fs);


if result == 0
    return
end

return

% Important input parameters:
% WedgeIndex, bendingAngle (thetaR-thetaS), minAngle (min(thetaS,
% wedgeIndex - thetaR)) or thetaS and (wedgeIndex - thetaR)

input = Geometry(step, shadowZone);

% Important input parameters:
% WedgeIndex, bendingAngle (thetaR-thetaS), minAngle (min(thetaS,
% wedgeIndex - thetaR)) or thetaS and (wedgeIndex - thetaR)

numInputs = length(input);
wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;

% Result template
rtemplate.ir = [];
rtemplate.tfmag = [];
rtemplate.tvec = [];
rtemplate.fvec = [];
rtemplate.tfcomplex = [];

result = repmat(rtemplate, 1, 1);

index1 = DataHash(input);
index2 = DataHash([wedgeLength, radiusS, radiusR, zS, zR, fs]);

index = [index1, index2];

% Create file info
mFile = mfilename('fullpath');
[inFilePath,fileStem] = fileparts(mFile);
fileStem = [fileStem, '_', num2str(index)];
savepath = ['results\', fileStem];
loadpath = ['results\', fileStem, '.mat'];

test = exist(loadpath, "file");

if test == 2
    load(loadpath, "result");
else
    for i = 1:numInputs
        wedgeIndex = input(i).wedge;
        thetaS = input(i).source + 0.01;
        thetaR = input(i).receiver - 0.01;
        if thetaR == 359.99
            thetaR = 359.98;
        end
        if wedgeIndex == 180
            wedgeIndex = 180.01;
        elseif wedgeIndex == 360
            wedgeIndex = 359.99;
        end
        [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec, result(i).tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs);
    end
end

save(savepath, "result");

% Input data
wedgeIndex = 345.01:10:355.01;
thetaW = 360 - wedgeIndex;
thetaS = 30;
thetaR = wedgeIndex - 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SingleWedge function

numResults = length(wedgeIndex);

[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex(1),thetaS,thetaR(1),radiusS,radiusR,zS,zR,fs);

tlen = length(tvec);
flen = length(fvec);

irall = zeros(numResults, tlen);
tfall = zeros(numResults, flen);

irall(1,:) = ir;
tfall(1,:) = tfmag;


% for i = 2:numResults
% 
%     [ir, tf, tvec, fvec] = SingleWedge(wedgeLength,wedgeIndex(i),thetaS,thetaR(i),radiusS,radiusR,zS,zR,fs);
% 
%     irall(i,:) = ir;
%     tfall(i,:) = tf;
% end

    T = 1 / fs;
    fc = 100;
    omega = 2 * pi * fc;
    K = omega * T;
    
    Hz.a0 = K;
    Hz.a1 = K;
    Hz.b0 = K + 2;
    Hz.b1 = K - 2;
    
    num = [Hz.a0 Hz.a1];
    den = [Hz.b0 Hz.b1];
    filter = dfilt.df1(num, den);

    [z, p, k] = tf2zpk(num, den);
    fc = 200;
    G = 2;
    Bpi = 10 ^ (G / 20);
    omega = 2 * pi * fc;
    K = omega * T;
    
    Hz.a0 = K + 2 * Bpi;
    Hz.a1 = K - 2 * Bpi;
    Hz.b0 = K + 2;
    Hz.b1 = K - 2;
    
    num = [Hz.a0 Hz.a1];
    den = [Hz.b0 Hz.b1];
    filter = dfilt.df1(num, den);

    [z, p, k] = tf2zpk(num, den);

%% For filters

% Start values
plpf = 0.998;
phsh = 0.5;
zlpf = -1;
zhsh = 0.8;
klpf = 0.002;
khsh = 30;

% in series
% pcomb = [plpf; phsh];
% zcomb = [zlpf; zhsh];
% kcomb = klpf * khsh;


% change k - overall level
% change p - higher the pole the higher the fc

z1 = -1;
z2= 0.5;
p1 = 0.5;
p2 = 0.5;
k1 = 0.1;

store = -1 + 2 * rand(5,1);
z1 = -1;
z1 = -1+0.5*abs(store(1));
z2= store(2);
p1 = store(3);
p2 = store(4);
k1 = abs(store(5)) / 10;

z = [z1; z2];
p = [p1; p2];
k = k1;

z1step = 0.1;
z2step = 0.1;
p1step = 0.1;
p2step = 0.1;
k1step = 0.01;
error = 100;
count = 0;
subcount = 0;

while error > 1 && count < 8
    improvement = false;
    [tfiir, ~] = IIRFilter(z, p, k, fs);
    error = Error(tfiir, tfmag, fvec);
    
    z2new = z2 + z2step;
    z(2,1) = z2new;
    [tfiir, ~] = IIRFilter(z, p, k, fs);
    errornew = Error(tfiir, tfmag, fvec);

    if errornew >= error
        z(2,1) = z2;
        z2step = -z2step;
    else
        improvement = true;
        z2 = z2new;
        error = errornew;
    end
    z1new = z1 + z1step;
    z(1,1) = z1new;
    [tfiir, ~] = IIRFilter(z, p, k, fs);
    errornew = Error(tfiir, tfmag, fvec);

    if errornew >= error
        z(1,1) = z1;
        z1step = -z1step;
    else
        improvement = true;
        z1 = z1new;
        error = errornew;
    end
    p2new = p2 + p2step;
    p(2,1) = p2new;
    [tfiir, fveciir] = IIRFilter(z, p, k, fs);
    errornew = Error(tfiir, tfmag, fvec);

    if errornew >= error
        p(2,1) = p2;
        p2step = -p2step;
    else
        improvement = true;
        p2 = p2new;
        error = errornew;
    end
    p1new = p1 + p1step;
    p(1,1) = p1new;
    [tfiir, fveciir] = IIRFilter(z, p, k, fs);
    errornew = Error(tfiir, tfmag, fvec);

    if errornew >= error
        p(1,1) = p1;
        p1step = -p1step;
    else
        improvement = true;
        p1 = p1new;
        error = errornew;
    end
    k1new = k1 + k1step;
    k = k1new;
    [tfiir, fveciir] = IIRFilter(z, p, k, fs);
    errornew = Error(tfiir, tfmag, fvec);

    if errornew >= error
        k = k1;
        k1step = -k1step;
    else
        improvement = true;
        k = k1new;
    end
    figure(2)
    semilogx(fveciir, tfiir)
    hold on
    semilogx(fvec, tfmag)
    hold off
    xlim([20, 20000])
    ylim([-40 0])

    if ~improvement
        subcount = subcount + 1;
    end
    if subcount > 6
        z2step = z2step / 2;
        p1step = p1step / 2;
        p2step = p2step / 2;
        k1step = k1step / 2;
        count = count + 1;
        subcount = 0;
    end
end

[b, a] = zp2tf(z, p, k);



% in parallel
% Not obvious how the poles and zeros change
% kcomb = klpf + khsh;

    % Not working as expected - works fine with matlab made filter
    % responses
%     fvecnorm = (fvec * 2 * pi) / fs;
%     [b, a] = invfreqz(tfcomplex, fvecnorm, 2, 2);
% 
%     [htest, f] = freqz(b, a, 2048);






H = zpk(z, p, k);


z2 = cell2mat(H.Z);
p2 = cell2mat(H.P);
k2 = H.K;
[b, a] = zp2tf(z, p, k);

H = dfilt.df1(b, a);

tfmag = filter(H, 1:100);

plot(tfmag, 1:length(tfmag))

testValues = zeros(flen, 1);

error = Error(testValues, tfmag, fvec);

figure(1)
plot(tvec, irall)

figure(2)
semilogx(fvec, tfall)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create file info
mFile = mfilename('fullpath');
[inFilePath,fileStem] = fileparts(mFile);

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
figure(3)
eddatafile = [inFilePath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
EDplotmodel(eddatafile,3,'figurewindow',3)

% Plot the irs
figure(4)
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
figure(5)
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


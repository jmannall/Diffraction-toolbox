%
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YTEA101
% Project Title: Particle Swarm Optimization Video Tutorial
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer and Instructor: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultLineMarkerSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SingleWedgeArray input data

wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;
wedge = 190:10:350;
bendingAngle = 190:10:350;
minAngle = 0:10:180;

[result, geometry] = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, wedge, bendingAngle, minAngle);

result = rmfield(result,'i');

%% Check for previous attempts

numResults = length(result);
count = 0;

%% Problem definition

problem.nVar = 5;       % Number of Unknown (Decision) Variables
problem.VarMin =  -0.99;  % Lower Bound of Decision Variables
problem.VarMax =  0.99;   % Upper Bound of Decision Variables

%% Parameters of PSO

params.MaxIt = 200;         % Maximum Number of Iterations
params.nPop = 200;          % Population Size (Swarm Size)
params.w = 0.5;             % Intertia Coefficient
params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Information

% Create file info
mFile = mfilename('fullpath');
[inFilePath,fileStem] = fileparts(mFile);
fileStem = [fileStem];

index = DataHash(result);
saveStem = [fileStem, 'Array_', num2str(index)];
savepath = ['results\', saveStem];
loadpath = ['results\', saveStem, '.mat'];

test = exist(loadpath, "file");

ptemplate.BestSol = [];
ptemplate.BestCosts = [];
ptemplate.tfmagiir = [];
ptemplate.fveciir = [];
ptemplate.z = [];
ptemplate.p = [];
ptemplate.k = [];
ptemplate.i = [];

pso = repmat(ptemplate, 1, 1);

saveThreshold = 2;
rerunThreshold = 0.5;

if test == 2
    load(loadpath, "pso");
    count = numResults;
    disp(['Result loaded from save']);
else
    for i = 1:numResults

        indexi = DataHash(result(i));
        oldCost = 10e999;   % To prevent count being added twice for new runs below save threshold
        rerun = true;
        
        % Create file info
        saveStemi = [fileStem, '_', num2str(indexi)];
        savepathi = ['results\psoIIRFilterConstrained\', saveStemi];
        loadpathi = ['results\psoIIRFilterConstrained\', saveStemi, '.mat'];
        
        test = exist(loadpathi, "file");

        if test == 2
            rerun = false;
            load(loadpathi, "BestSol", "BestCosts");
            oldCost = BestCosts(end);
            oldBestSol = BestSol;
            oldBestCosts = BestCosts;
            oldk = BestSol.Position(5);
            if oldCost > rerunThreshold
                rerun = true;
                saveThreshold = oldCost;
            end
        end
        if rerun
            %% Problem Definiton - cost function
            
            problem.CostFunction = @(x) FilterError(x, result(i).tfmag, fs);  % Cost Function
            problem.ConstraintFunction = @(x) FilterConstraints(x); % Constraint Function
            
            %% Calling PSO
            
            out = PSOWithIIRConstraints(problem, params);
            
            BestSol = out.BestSol;
            BestCosts = out.BestCosts;
        end
        
        %% Results
        
        [filterStable, tfmagiir, fveciir, b, a] = psoProcessResults(BestSol, fs);
        
        [z, p, k] = tf2zpk(b, a);
    
        if filterStable
            if BestCosts(end) <= saveThreshold
                save(savepathi, "BestSol", "BestCosts");
                pso(i).BestSol = BestSol;
                pso(i).BestCosts = BestCosts;
                pso(i).tfmagiir = tfmagiir;
                pso(i).fveciir = fveciir;
                pso(i).z = z;
                pso(i).p = p;
                pso(i).k = k;
                count = count + 1;
            end
            if oldCost < BestCosts(end)
                count = count + 1;
            end
        end

%         if filterStable
%             if BestCosts(end) <= oldCost || oldk < 0
%                 save(savepathi, "BestSol", "BestCosts");
%                 pso(i).BestSol = BestSol;
%                 pso(i).BestCosts = BestCosts;
%                 pso(i).tfmagiir = tfmagiir;
%                 pso(i).fveciir = fveciir;
%                 pso(i).z = z;
%                 pso(i).p = p;
%                 pso(i).k = k;
%             end
%             if BestCosts(end) < threshold
%                 count = count + 1;
%             end
%         end


        
        f1 = figure(1);
        movegui(f1,'northwest');
        % plot(BestCosts, 'LineWidth', 2);
        semilogy(BestCosts);
        xlabel('Iteration');
        ylabel('Best Cost');
        grid on;
        
        f2 = figure(2);
        movegui(f2,'north');
        semilogx(fveciir, tfmagiir)
        hold on
        semilogx(result(i).fvec, result(i).tfmag)
        hold off
        xlim([20, 20000])
        ylim([-40 0])
        
        ts = 1 / fs;
        sys = tf(b, a, ts);

        f3 = figure(3);
        movegui(f3,'northeast')
        h = pzplot(sys);
        grid on
    
    end
end

%% Save result

disp([num2str(count) ' out of ' num2str(numResults) ' filters found within threshold']);

if count == numResults
    save(savepath, "pso");
end

%% Process data

x = [pso.BestSol];
resultCost = [x.Cost];
resultPosition = reshape([x.Position], problem.nVar, numResults);
worstCost = max(resultCost);
averageCost = sum(resultCost) / numResults;

z = resultPosition([1,3],:);
p = resultPosition([2,4],:);
k = resultPosition(5,:);

[wedge, ~, iw] = unique(geometry.wedge);
[bendingAngle, ~, ibA] = unique(geometry.bendingAngle);
[minAngle, ~, imA] = unique(geometry.minAngle);

numWedges = length(wedge);
numBendingAngles = length(bendingAngle);
numMinAngles = length(minAngle);

close all

wbA = [geometry.wedge, geometry.bendingAngle];
[WBA, ~, iwbA] = unique(wbA, 'rows');
wmA = [geometry.wedge, geometry.minAngle];
[WMA, ~, iwmA] = unique(wmA, 'rows');
bAmA = [geometry.bendingAngle, geometry.minAngle];
[BAMA, ~, ibAmA] = unique(bAmA, 'rows');

%% Define filters used

zAll = [resultPosition(1,:); resultPosition(3,:)];
pAll = [resultPosition(2,:); resultPosition(4,:)];
kAll = [resultPosition(5,:)];
tfmagiirAll = [pso.tfmagiir];

target = tfmagiirAll(1,:) - 3;
y = (10.^(target / 10) ./ kAll.^2);
a = (y .* (pAll(2,:).^2 + 1) - zAll(1,:).^2 - 1) ./ (2 .* (pAll(2,:) .* y - zAll(1,:)));
b = sqrt(1 - a.^2);
zNum = a+b*1i;
omegac = angle(zNum);
fclpf = omegac / (2 * pi) * fs;

target = tfmagiirAll(end / 2,:) - 3;
y = (10.^(target / 10) ./ kAll.^2);
a = (y .* (pAll(2,:).^2 + 1) - zAll(1,:).^2 - 1) ./ (2 .* (pAll(2,:) .* y - zAll(1,:)));
b = sqrt(1 - a.^2);
zNum = a+b*1i;
omegac = angle(zNum);
fchsh = omegac / (2 * pi) * fs;

B0 = (1 - pAll(1,:)) .* (1 + zAll(2,:)) ./ ((1 + pAll(1,:)) .* (1 - zAll(2,:)));

disp('section complete');

%% Regression

w = geometry.wedge;
bA = geometry.bendingAngle;
mA = geometry.minAngle;
w2 = w.^2;
bA2 = bA.^2;
mA2 = mA.^2;
%x = [w, bA, mA, w .* bA, w .* mA, bA .* mA, w .* bA .* mA, w2, bA2, mA2, w2 .* bA, w2 .* mA, bA2 .* w, bA2 .* mA, mA2 .* w, mA2 .* bA, w2 .* bA2, w2 .* mA2, bA2 .* mA2, w2 .*mA2.*bA2];

%x = [w, bA, mA, w.*mA, w.*bA, bA.*mA, mA.^2, bA.^2, w.^2, bA ./ w, mA ./ w, bA ./ (mA + 1), tand(bA+1), tand(mA+1), tand(w+1), cosd(bA), cosd(mA), sind(w), tanh(180 ./ (2 * pi * bA)), cosd(w), bA ./ (w - mA - 180)];
%x = [w, bA, mA, w.*mA, w.*bA, bA.*mA, mA.^2, bA.^2, w.^2, tanh(1 ./ (2 * pi * (mA + 1))), cosd(mA), tanh(180 ./ (2 * pi * w)), tanh(180 ./ (2 * pi * bA)), bA ./ (w - mA - 180)];
%x = [w, bA, mA, w.*mA, w.*bA, bA.*mA, tanh(1 ./ (2 * pi * (mA + 1))), tanh(180 ./ (2 * pi * w)), tanh(180 ./ (2 * pi * bA)), bA ./ (w - mA - 180)];
%x = [w, bA, mA];
%x = [sind(w), sind(bA), cosd(mA)];
%x = [w, bA, mA, bA.^2, bA .* mA, w .* mA, w .* mA .* bA, sind(w) .* cos(mA), sind(w) .* sind(bA), sind(bA) .* cosd(mA), sind(w), sind(bA), cosd(mA), tanh(1 ./ (2 * pi * mA + 0.2)), tanh(180 ./ (2 * pi * w)), tanh(180 ./ (2 * pi * bA)), bA ./ (w - mA - 180)];
%x = [tanh(180 ./ (2 * pi * bA)), tanh(180 ./ (2 * pi * w)), sind(bA), sind(bA) .* cosd(mA), sind(w) .* sind(bA), cosd(mA), bA, sind(w), mA, bA ./ (w - mA - 180), w, bA.^2, bA .* mA, w .* mA, w .* mA .* bA];
%x = [w, bA, mA, bA.^2, w .* bA, w .* mA, bA .* mA, w .* bA .* mA, sind(w), sind(bA), cosd(mA), tanh(1 ./ (2 * pi * mA + 0.2)), tanh(180 ./ (2 * pi * w)), tanh(180 ./ (2 * pi * bA)), bA ./ (w - mA - 180)];
x = [w, bA, mA, bA.^2, w ./ bA, bA .* mA, w .* mA, mA .^ 2, w .* bA, sqrt(abs(sind(bA))), tanh(180 ./ (2 * pi * w)), tanh(180 ./ (2 * pi * bA)), bA ./ (w - mA - 180), w .* mA .* bA, sind(bA)];

% w = w * pi / 180;
% bA = bA * pi / 180;
% mA = mA * pi / 180;
% x = [w, bA, mA, bA.^2, w ./ bA, bA .* mA, w .* mA, mA .^ 2, w .* bA, sqrt(abs(sin(bA))), tanh(1 ./ (2 * w)), tanh(1 ./ (2 * bA)), bA ./ (w - mA - pi), w .* mA .* bA, sin(bA)];


input = transpose(resultPosition);
z1 = input(:,1);
p1 = input(:,2);
z2 = input(:,3);
p2 = input(:,4);
k = input(:,5);

X = [ones(numResults, 1), x];
Y = input;
disp('break');
[beta,Sigma,E,CovB,logL] = mvregress(X,Y);
%[beta,Sigma,E,CovB,logL] = mvregress(X,Y, 'outputfcn', @(x, y, z) Test(x, y, z, X, [result.tfmag], fs));
averageError = sum(E.^ 2,1) / length(E);
worstError = max(E,[],1);

PosE = sum(E,2);
averagePosError = sum(PosE.^ 2,1)  / length(PosE);
worstPosError = max(PosE,[],1);

B = beta
sig = sum(abs(B),2);
prediction = ([ones(numResults,1),x]*B)';

zPrediction = [prediction(1,:); prediction(3,:)];
pPrediction = [prediction(2,:); prediction(4,:)];
kPrediction = [prediction(5,:)];

MSE = ZPKError(B, X, [result.tfmag], fs);

disp('Section complete');

%% Compare filters
close all

for i = 1:numResults
        
        input = zPrediction(i);
        [tfmagPrediction, fvecPrediction] = IIRFilter(zPrediction(:,i), pPrediction(:,i), kPrediction(i), fs);
        
        f2 = figure;
        movegui(f2, 'northeast')
        semilogx(pso(i).fveciir, pso(i).tfmagiir)
        hold on
        semilogx(result(i).fvec, result(i).tfmag)
        semilogx(fvecPrediction, tfmagPrediction)
        hold off
        xlim([20, 20000])
        ylim([-40 0])
        
        saveas(gcf,['figures/filter_predictions/', 'w_', num2str(geometry.wedge(i)), '_bA_', num2str(geometry.bendingAngle(i)), '_mA_', num2str(geometry.minAngle(i)), '.png'])
        close all
end

%% Plot graphs

% Change minAngle
plotFilterParameters(WBA, iwbA, geometry.minAngle, fclpf, fchsh, B0, 'w', 'bA');
% Change bendingAngle
plotFilterParameters(WMA, iwmA, geometry.bendingAngle, fclpf, fchsh, B0, 'w', 'mA');
% Change wedge
plotFilterParameters(BAMA, ibAmA, geometry.wedge, fclpf, fchsh, B0, 'bA', 'mA');

%% Plot graphs

% Change minAngle
plotFilterPZK(WBA, iwbA, geometry.minAngle, z, p, k, 'w', 'bA', zPrediction, pPrediction, kPrediction);
% Change bendingAngle
plotFilterPZK(WMA, iwmA, geometry.bendingAngle, z, p, k, 'w', 'mA', zPrediction, pPrediction, kPrediction);
% Change wedge
plotFilterPZK(BAMA, ibAmA, geometry.wedge, z, p, k, 'bA', 'mA', zPrediction, pPrediction, kPrediction);

return

%% OLD code
p = [pso.p];
z = [pso.z];
k = [pso.k];
pole1 = p(1,:);
pole2 = p(2,:);
zero1 = z(1,:);
zero2 = z(2,:);
gain = k;

[~, ~, Iw] = unique(w);
[~, ~, ImA] = unique(mA);
[~, ~, IbA] = unique(bA);

[pdata1, pdata2, zdata1, zdata2, kdata] = deal(-5 * ones(size(w)));

value = [iw, imA, ibA];
match = [Iw, ImA, IbA];

[q, idx] = ismember(value, match, 'rows');

for i = 1:length(idx)
    pdata1(idx(i)) = pole1(i);
    pdata2(idx(i)) = pole2(i);
    zdata1(idx(i)) = zero1(i);
    zdata2(idx(i)) = zero2(i);
    kdata(idx(i)) = gain(i);
end

test1 = reshape(squeeze(mA(:,19,:)),1,[]);
test2 = reshape(squeeze(bA(:,19,:)),1,[]);
test3 = reshape(squeeze(kdata(:,19,:)),1,[]);

index = find(test3 == -5);

test1(index) = [];
test2(index) = [];
test3(index) = [];

%% functions

% syms wsym mAsym bAsym
% 
% func_k = a*exp(b*bAsym) + c*exp(d*bAsym) + e*mAsym + f*bAsym + g

disp('Functions');

index1 = find(kdata ~= -5);
index2 = find(bA ~= 180);

[q, idx] = ismember(index1, index2);

count = 1;
index = [];
for i = 1:length(idx)
    if q(i)
        index(count,1) = index1(i);
        count = count + 1;
    end
end

values = zeros(length(minAngle),length(bendingAngle),length(wedge),3);
values(:,:,:,1) = bA;
values(:,:,:,2) = mA;
values(:,:,:,3) = w;
num = numel(bA);

input1 = values(index);
input2 = values(num + index);
input3 = values(2*num + index);

ydatak = reshape(kdata(index),[],1);
ydatap1 = reshape(pdata1(index),[],1);
ydatap2 = reshape(pdata2(index),[],1);
ydataz1 = reshape(zdata1(index),[],1);
ydataz2 = reshape(zdata2(index),[],1);
xdata = zeros(size(ydatak));

kfunc = @(x,xdata)KFunction(x,xdata,input1,input2,input3);
p1func = @(x,xdata)P1Function(x,xdata,input1,input2,input3);
p2func = @(x,xdata)P2Function(x,xdata,input1,input2,input3);
z1func = @(x,xdata)Z1Function(x,xdata,input1,input2,input3);
z2func = @(x,xdata)Z2Function(x,xdata,input1,input2,input3);


x0k = [0.1,-0.4,0,-7,0,0,0];
[kparams, kresnorm] = lsqcurvefit(kfunc,x0k,xdata,ydatak);
x0 = [1,1,1,0];
x0p1 = [0.1,-0.4,0,-7,0,0,0];
x0p2 = [0.1,-0.4,0,-7,0,0,0];
x0z2 = [0.1,-0.4,0,-7,0,0,0,1,1,1,1,1,1,1];
[p1params, p1resnorm] = lsqcurvefit(p1func,x0p1,xdata,ydatap1);
[p2params, p2resnorm] = lsqcurvefit(p2func,x0p2,xdata,ydatap2);
[z1params, z1resnorm] = lsqcurvefit(z1func,x0p1,xdata,ydataz1);
[z2params, z2resnorm] = lsqcurvefit(z2func,x0z2,xdata,ydataz2);

close all

times = 1:1:length(xdata);
figure(1)
plot(times,ydatak,'ko',times,kfunc(kparams,xdata),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')
figure(2)
plot(times,ydatak-kfunc(kparams,xdata))
grid on

figure(3)
plot(times,ydatap1,'ko',times,p1func(p1params,xdata),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')
figure(4)
plot(times,ydatap1-p1func(p1params,xdata))
grid on

figure(5)
plot(times,ydatap2,'ko',times,p2func(p2params,xdata),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')
figure(6)
plot(times,ydatap2-p2func(p2params,xdata))
grid on

figure(7)
plot(times,ydataz1,'ko',times,z1func(z1params,xdata),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')
figure(8)
plot(times,ydataz1-z1func(z1params,xdata))
grid on

figure(9)
plot(times,ydataz2,'ko',times,z2func(z2params,xdata),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')
figure(10)
plot(times,ydataz2-z2func(z2params,xdata))
grid on

%% Check result
close all

for i = 200:numResults
        
    iircoeff = [pso(i).BestSol];
    position = [iircoeff.Position];
    pole = [position(1); position(2)];
    zero = [position(3); position(4)];
    gain = [position(5)];



    bA = geometry(i).bendingAngle;
    mA = geometry(i).minAngle;
    w = geometry(i).wedge;

    %p1 = P1Function(p1params,[0,0,0],bA,mA,w);
    p1 = 0.95;
    p2 = P2Function(p2params,[0,0,0],bA,mA,w);
    %z1 = Z1Function(z1params,[0,0,0],bA,mA,w);
    z1 = 0.95;
    z2 = Z2Function(z2params,[0,0,0],bA,mA,w);
    k = KFunction(kparams,[0,0,0],bA,mA,w);
    
    polenew = [p1;p2];
    zeronew = [z1;z2];
    gainnew = k;

    [tfmagiir fveciir] = IIRFilter(zeronew,pole,gainnew,fs);
    [tfmagiir2 fveciir2] = IIRFilter(zero,pole,gain,fs);

    f2 = figure(2);
    movegui(f2,'north');
    semilogx(fveciir2, tfmagiir2)
    hold on
    semilogx(result(i).fvec, result(i).tfmag)
    semilogx(fveciir, tfmagiir)
    hold off
    xlim([20, 20000])
    ylim([-40 0])

    [a,b] = zp2tf(pole, zero, gain);
    [anew,bnew] = zp2tf(pole, zeronew, gainnew);
    [poletest, zerotest, gaintest] = tf2zpk(anew,bnew);
    ts = 1/fs;
    sys = tf(a,b,ts);
    sysnew = tf(anew,bnew,ts);
        
    f3 = figure(3);
    movegui(f3,'northeast')
    h = pzplot(sys);
    grid on
    f4 = figure(4);
    movegui(f4,'southeast')
    g = pzplot(sysnew);
    grid on
end


for i = 200:numResults
    bA = geometry(i).bendingAngle;
    mA = geometry(i).minAngle;
    w = geometry(i).wedge;

    %p1 = P1Function(p1params,[0,0,0],bA,mA,w);
    p1 = 0.95;
    p2 = P2Function(p2params,[0,0,0],bA,mA,w);
    %z1 = Z1Function(z1params,[0,0,0],bA,mA,w);
    z1 = 0.95;
    z2 = Z2Function(z2params,[0,0,0],bA,mA,w);
    k = KFunction(kparams,[0,0,0],bA,mA,w);
    
    [tfmagiir fveciir] = IIRFilter([z1;z2],[p1;p2],k,fs);

    f2 = figure(2);
    movegui(f2,'north');
    semilogx(fveciir, tfmagiir)
    hold on
    semilogx(result(i).fvec, result(i).tfmag)
    hold off
    xlim([20, 20000])
    ylim([-40 0])
end

close all

%tiledlayout(10,19) % Change to plot single variable change for other two variables kept constant.
for i = 1:length(wedge)
    figure('Position',[10 250 1900 560])
    tiledlayout(1,3)
    nexttile
    plot3(squeeze(mA(:,i,:)),squeeze(bA(:,i,:)),squeeze(pdata1(:,i,:)));
    hold on
    plot3(squeeze(mA(:,i,:)),squeeze(bA(:,i,:)),squeeze(pdata2(:,i,:)));
    hold off
    zlim([-1, 1]);
    title('Poles');
    nexttile
    plot3(squeeze(mA(:,i,:)),squeeze(bA(:,i,:)),squeeze(zdata1(:,i,:)));
    hold on
    plot3(squeeze(mA(:,i,:)),squeeze(bA(:,i,:)),squeeze(zdata2(:,i,:)));
    hold off
    zlim([-1, 1]);
    title('Zeros');
    nexttile
    plot3(squeeze(mA(:,i,:)),squeeze(bA(:,i,:)),squeeze(kdata(:,i,:)));
    zlim([0, 1]);
    title('Gain');
end

sliceomatic(data,w,mA,bA)

wedge = [geometry.wedge];
source = [geometry.source];
receiver = [geometry.receiver];
bendingAngle = [geometry.bendingAngle];
minAngle = [geometry.minAngle];

[~, ~, icw] = unique(wedge);
[~, ~, ics] = unique(source);
[~, ~, icr] = unique(receiver);
[~, ~, icb] = unique(bendingAngle);
[~, ~, icm] = unique(minAngle);

numWedges = max(icw);
numSources = max(ics);
numReceivers = max(icr);
numBendingAngles = max(icb);
numMinAngles = max(icm);

itemplate.i = [];
itemplate.w = [];
itemplate.s = [];
itemplate.r = [];

ialltemplate.wedge = [];
ialltemplate.source = [];
ialltemplate.receiver = [];

indexw = repmat(itemplate, 1, 1);
indexs = repmat(itemplate, 1, 1);
indexr = repmat(itemplate, 1, 1);
indexall = repmat(itemplate, 1, 1);

for i = 1:numWedges
    indexw(i).i = find(icw == i);
end
for i = 1:numSources
    indexs(i).i = find(ics == i);
end
for i = 1:numReceivers
    indexr(i).i = find(icr == i);
end
for i = 1:numBendingAngles
    indexb(i).i = find(icb == i);
end
for i = 1:numMinAngles
    indexm(i).i = find(icm == i);
end


for i = 1:numWedges
    for j = 1:numSources
        for k = 1:numReceivers
        indexall(numWedges + 1,j,k).i = intersect(indexs(j).i, indexr(k).i);
        indexall(i,numSources + 1,k).i = intersect(indexw(i).i, indexr(k).i);
        indexall(i,j,numReceivers + 1).i = intersect(indexw(i).i, indexs(j).i);
        [indexall(numWedges + 1,j,k).w, indexall(i,numSources + 1,k).s, indexall(i,j,numReceivers + 1).r] = deal("variable");
        [indexall(i,j,numReceivers + 1).w, indexall(i,numSources + 1,k).w] = deal(wedge(indexw(i).i(1)));
        [indexall(numWedges + 1,j,k).s, indexall(i,j,numReceivers + 1).s] = deal(source(indexs(j).i(1)));
        [indexall(numWedges + 1,j,k).r, indexall(i,numSources + 1,k).r] = deal(receiver(indexr(k).i(1)));
        end
    end
end

indexall = [];

for i = 1:numWedges
    for j = 1:numBendingAngles
        for k = 1:numMinAngles
        indexall(numWedges + 1,j,k).i = intersect(indexb(j).i, indexm(k).i);
        indexall(i,numBendingAngles + 1,k).i = intersect(indexw(i).i, indexm(k).i);
        indexall(i,j,numMinAngles + 1).i = intersect(indexw(i).i, indexb(j).i);
        [indexall(numWedges + 1,j,k).w, indexall(i,numBendingAngles + 1,k).s, indexall(i,j,numMinAngles + 1).r] = deal("variable");
        [indexall(i,j,numMinAngles + 1).w, indexall(i,numBendingAngles + 1,k).w] = deal(wedge(indexb(i).i(1)));
        [indexall(numWedges + 1,j,k).s, indexall(i,j,numMinAngles + 1).s] = deal(source(indexb(j).i(1)));
        [indexall(numWedges + 1,j,k).r, indexall(i,numBendingAngles + 1,k).r] = deal(receiver(indexm(k).i(1)));
        end
    end
end

indexing = squeeze(indexall(:,numBendingAngles + 1,:));

variable = "bendingangle";

p = [pso.p];
p1 = p(1,:);
p2 = p(2,:);

[newPoles1, newPoles2] = deal(5 * ones(numBendingAngles, numMinAngles));
for i = 1:numResults
    newPoles1(icb(i),icm(i)) = p(1,i);
    newPoles2(icb(i),icm(i)) = p(2,i);
end

a = unique(minAngle);
b = unique(bendingAngle);

figure(7)
mesh(a,b, newPoles1)

figure(8)
mesh(a,b, newPoles2)

figure(9)
plot3(minAngle, bendingAngle, p(1,:),'x')
hold on
plot3(minAngle, bendingAngle, p(2,:),'o')

PlotFigures(pso, bendingAngle, indexing, variable, geometry);

figure(5)
scatter([geometry.source], scale)

figure(6)
scatter([geometry.receiver], scale)

figure(7)
scatter([geometry.bendingAngle], scale)

figure(8)
scatter([geometry.minAngle], scale)

function output = KFunction(p,xdata,x,y,z)
    mat = p(1).*exp(p(2).*x)+p(3).*exp(p(4).*x)+p(5).*y+p(6).*z+p(7);
    output = reshape(mat,[numel(x),1]);
end

function output = P1Function(p,xdata,x,y,z)
    mat = p(1).*x+p(2).*y+p(3).*z+p(4).*x.^2+p(5).*y.^2+p(6).*z.^2+p(7);
    output = reshape(mat,[numel(x),1]);
end

function output = P2Function(p,xdata,x,y,z)
    mat = p(1).*x+p(2).*y+p(3).*z+p(4).*x.^2+p(5).*y.^2+p(6).*z.^2+p(7);
    output = reshape(mat,[numel(x),1]);
end

function output = Z1Function(p,xdata,x,y,z)
    mat = p(1).*x+p(2).*y+p(3).*z+p(4).*x.^2+p(5).*y.^2+p(6).*z.^2+p(7);
    output = reshape(mat,[numel(x),1]);
end

function output = Z2Function(p,xdata,x,y,z)
    mat = (p(1).*x+p(2).*y+p(3).*z+p(4).*x.^2+p(5).*y.^2+p(6).*z.^2+p(7))./(p(8).*x+p(9).*y+p(10).*z+p(11).*x.^2+p(12).*y.^2+p(13).*z.^2+p(14));
    output = reshape(mat,[numel(x),1]);
end

function Objective = Test(current, input, text, X, tfvalue, fs)
    disp(['Iteration: ', num2str(input.iteration)]);
    current = reshape(current, 16, 5);
    MSE = ZPKError(current, X, tfvalue, fs);
    disp(['MSE: ', num2str(MSE)])
    Objective = true;
    if strcmp(text,'done')
        disp('Iterations complete');
        Objective = false;
    end
end

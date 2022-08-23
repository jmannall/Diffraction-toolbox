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

% Input data
wedgeLength = 20;
wedgeIndex = 275.01;
thetaW = 360 - wedgeIndex;
thetaS = 10.01;
thetaR = 250.01;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;

fs = 96000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SingleWedge function

numResults = length(wedgeIndex);

% [ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,fs);

% SingleWedgeArray function

wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;

step = 20;
minw = 180;
maxw = 360;
shadowZone = false;

[result, geometry] = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, step, shadowZone, minw, maxw);

%% Check for previous attempts

numResults = length(result);
count = 0;

%% Problem definition

problem.nVar = 5;       % Number of Unknown (Decision) Variables
problem.VarMin =  -1;  % Lower Bound of Decision Variables
problem.VarMax =  1;   % Upper Bound of Decision Variables
problem.split = 5;      % 3 if FilterErrorConj / 5 if FilterError

%% Parameters of PSO

params.MaxIt = 100;         % Maximum Number of Iterations
params.nPop = 100;          % Population Size (Swarm Size)
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

if test == 2
    load(loadpath, "pso");
    count = numResults;
    disp(['Result loaded from save']);
else
    for i = 1:numResults
        indexi = DataHash(result(i));
        
        % Create file info
        saveStemi = [fileStem, '_', num2str(indexi)];
        savepathi = ['results\', saveStemi];
        loadpathi = ['results\', saveStemi, '.mat'];
        
        test = exist(loadpathi, "file");

        if test == 2
            load(loadpathi, "BestSol", "BestCosts");
        else
            %% Problem Definiton - cost function
            
            problem.CostFunction = @(x) ParameterError(x, result(i).tfmag, geometry(i), fs);  % Cost Function
            
            %% Calling PSO
            
            out = PSO(problem, params);
            
            BestSol = out.BestSol;
            BestCosts = out.BestCosts;
        end
        
        %% Results
        
        [filterStable, tfmagiir, fveciir, b, a] = psoProcessResults(BestSol, fs);
        
        costEnd = length(BestCosts);
        [z, p, k] = tf2zpk(b, a);
    
        if filterStable && BestCosts(costEnd) < 1
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

disp([num2str(count) ' out of ' num2str(numResults) ' filters found within threshold']);

if count == numResults
    save(savepath, "pso");
end

return

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

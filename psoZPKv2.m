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

disp('Create training data');

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

problem.VarMin =  -1;  % Lower Bound of Decision Variables
problem.VarMax =  1;   % Upper Bound of Decision Variables

%% Parameters of PSO

params.MaxIt = 40;         % Maximum Number of Iterations
params.nPop = 200;          % Population Size (Swarm Size)
params.w = 0.5;             % Intertia Coefficient
params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Information

% Create file info
mFile = mfilename('fullpath');
[inFilePath,fileStem] = fileparts(mFile);

modelNum = 6;

index = DataHash(result);
saveStem = [fileStem, '_model', num2str(modelNum)];
savepath = ['results\', saveStem];
loadpath = ['results\', saveStem, '.mat'];

w = geometry.wedge;
bA = geometry.bendingAngle;
mA = geometry.minAngle;

input = ([ones(numResults, 1), sind(w / 2), sind(bA / 2), cosd(mA)])';

startingValues = ([0.99, 0.9, 0.5, -0.5, 0.2;
    -0.2, 0, 0, -0.2, 0;
    -0.2, -0.2, -0.1, 0.5, 0.1;
    -0.8, 0.1, -0.5, 0, -0.3])';
% startingValues = [0.990000000000000	0.0151775062996523	-0.399390163720326	-0.129178938122147
% 0.685078624238558	0.0111212824039604	-0.0830819731343385	0.321193163210682
% 0.520100350459114	-0.537769567431370	-0.285561796943554	-1.17780309293581
% -0.872849183466677	-0.00748942501053247	0.717877023069611	-1.11072074945777
% 0.0434710692424342	0.00181621204324149	0.0518723307746558	-0.0389711453278028];

problem.nVar = [size(startingValues)];       % Number of Unknown (Decision) Variables

test = exist(loadpath, "file");

if test == 2
    load(loadpath, "BestSol", "BestCosts");
    count = numResults;
    disp('Result loaded from save');
else
    disp('Begin PSO');
    problem.CostFunction = @(x) ZPKv2Error(x, input, [result.tfmag], fs);  % Cost Function
    problem.ConstraintFunction = @(x) ZPKv2Constraints(x, input); % Constraint Function
    problem.InitialiseFunction = @(x) ZPKv2Initialise(x, startingValues); % Initialisation Function
    
    %% Calling PSO
    
    out = PSOWithConstraints(problem, params);
    
    BestSol = out.BestSol;
    BestCosts = out.BestCosts;

    save(savepath, "BestSol", "BestCosts", "input")
end

%% Plot data

f1 = figure(1);
movegui(f1,'southwest');
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

prediction = BestSol.Position * input;

z = [prediction(1,:); prediction(3,:)];
p = [prediction(2,:); prediction(4,:)];
k = [prediction(5,:)];

for i = 1:numResults
    [tfmagModel, fvecModel] = IIRFilter(z(:,i), p(:,i), k(i), fs);

    f2 = figure(2);
    movegui(f2,'south');
    semilogx(fvecModel, tfmagModel)
    hold on
    semilogx(result(i).fvec, result(i).tfmag)
    hold off
    xlim([20, 20000])
    ylim([-40 0])
        
    [b, a] = zp2tf(z(:,i), p(:,i), k(i));
    ts = 1 / fs;
    sys = tf(b, a, ts);

    f3 = figure(3);
    movegui(f3,'southeast')
    h = pzplot(sys);
    grid on
end
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
set(groot, 'defaultLineMarkerSize', 20)

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
shadowZone = true;

result = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, step, shadowZone);

%% Check for previous attempts

numResults = length(result);
count = 0;

%% Problem definition

problem.nVar = 5;       % Number of Unknown (Decision) Variables - Return to 5. All filters end up with 0 imaginary anyway?
problem.VarMin =  -1;  % Lower Bound of Decision Variables
problem.VarMax =  1;   % Upper Bound of Decision Variables

%% Parameters of PSO

params.MaxIt = 200;         % Maximum Number of Iterations
params.nPop = 400;          % Population Size (Swarm Size)
params.w = 0.5;             % Intertia Coefficient
params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Information

for i = 1:numResults
    index = DataHash(result(i));
    
    % Create file info
    mFile = mfilename('fullpath');
    [inFilePath,fileStem] = fileparts(mFile);
    fileStem = [fileStem, '_', num2str(index)];
    savepath = ['results\', fileStem];
    loadpath = ['results\', fileStem, '.mat'];
    
    test = exist(loadpath, "file");
    
    if test == 2
        load(loadpath, "BestSol", "BestCosts");
    else
        %% Problem Definiton - cost function
        
        problem.CostFunction = @(x) FilterError(x, result(i).tfmag, fs);  % Cost Function
        
        %% Calling PSO
        
        out = PSO(problem, params);
        
        BestSol = out.BestSol;
        BestCosts = out.BestCosts;
    end
    
    %% Results
    
    [filterStable, tfiir, fveciir, b, a] = psoProcessResults(BestSol, fs);
    
    costEnd = length(BestCosts);

    if filterStable && BestCosts(costEnd) < 0.885
        save(savepath, "BestSol", "BestCosts");
        count = count + 1;
    end
    
    figure(1)
    % plot(BestCosts, 'LineWidth', 2);
    semilogy(BestCosts);
    xlabel('Iteration');
    ylabel('Best Cost');
    grid on;
    
    figure(2)
    semilogx(fveciir, tfiir)
    hold on
    semilogx(result(i).fvec, result(i).tfmag)
    hold off
    xlim([20, 20000])
    ylim([-40 0])
    
    ts = 1 / fs;
    sys = tf(b, a, ts);
    
    figure(3)
    h = pzplot(sys);
    grid on

end

disp([num2str(count) ' out of ' num2str(numResults) ' filters found within threshold']);


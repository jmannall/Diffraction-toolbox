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
wedgeIndex = 275.01:10:355.01;
thetaW = 360 - wedgeIndex;
thetaS = 10;
thetaR = wedgeIndex - 40;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
wedgeSize = max(radiusS, radiusR);

fs = 96000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SingleWedge function

numResults = length(wedgeIndex);

[ir, tfvalue, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength,wedgeIndex(1),thetaS,thetaR(1),radiusS,radiusR,zS,zR,fs);

%% Problem Definiton

problem.CostFunction = @(x) FilterError(x, tfvalue, fs);  % Cost Function
problem.nVar = 5;       % Number of Unknown (Decision) Variables - Return to 5. All filters end up with 0 imaginary anyway?
problem.VarMin =  -1;  % Lower Bound of Decision Variables
problem.VarMax =  1;   % Upper Bound of Decision Variables

%% Parameters of PSO

params.MaxIt = 50;        % Maximum Number of Iterations
params.nPop = 100;           % Population Size (Swarm Size)
params.w = 0.5;               % Intertia Coefficient
params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
params.c1 = 2;              % Personal Acceleration Coefficient
params.c2 = 2;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

%% Calling PSO

out = PSO(problem, params);

BestSol = out.BestSol;
BestCosts = out.BestCosts;

%% Results

z = [BestSol.Position(1);BestSol.Position(2)];
p = [BestSol.Position(3);BestSol.Position(4)];
k = BestSol.Position(5);

[b, a] = zp2tf(z, p, k);
[tfiir, fveciir] = IIRFilter(z, p, k, fs);

stabilityCheck = abs(p);

if stabilityCheck <= 1
    filterStable = true;
    str = 'is';
else
    filterStable = false;
    str = 'is not';
end

disp(['Filter ' str ' stable.']);

figure(1)
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

figure(2)
semilogx(fveciir, tfiir)
hold on
semilogx(fvec, tfvalue)
hold off
xlim([20, 20000])
ylim([-40 0])

ts = 1 / fs;
sys = tf(b, a, ts);

figure(3)
h = pzplot(sys);
grid on




close all
clear all

fs = 48000;

wedgeIndex = 350;
bendingAngle = 60;
minAngle = 20;

S.w = wedgeIndex;
S.bA = bendingAngle;
S.mA = minAngle;

plotFcn = @(S) BTMInf(S);
sliderFcn = @(S) CreateBTMSliders(S);
updateFcn = @UpdateBTM;

CreateInteractivePlot(S, plotFcn, sliderFcn, updateFcn, fs);











close all
clear all

numPillars = 13;
gratingWidth = 0.5;
pillarHeight = 10;

[corners, planeCorners, planeRigid, source, receiver] = CreateDiffractionGratingGeometry(numPillars, gratingWidth, pillarHeight);

    
    
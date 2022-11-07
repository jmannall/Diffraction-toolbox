function [geometry, input] = RandomGeometryWedge(numObservations)

    % Geometry template
    gtemplate.wedgeIndex = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);

    epsilon = 1e-5;

    %% Wedge index
    numTriObservations = 0.5 * numObservations;
    numUniObservations = numObservations - numTriObservations;

    const = ones(numObservations, 1);

    wedgeIndex = [181, 360];
    wI = [RandomTriangularDistribution(wedgeIndex, false, numTriObservations); RandomUniformDistribution(wedgeIndex, numUniObservations)];
    
    wedgeIndex = [wI, (360 - epsilon) * const];
    wI = RandomTriangularDistribution(wedgeIndex, true, numObservations);
    
    %% Bending angle and minimum angle
    bendingAngle = [(180 + epsilon) * const, wI - epsilon];
    bA = RandomTriangularDistribution(bendingAngle, false, numObservations);
    
    minAngle = [epsilon * const, (wI - bA) / 2];
    mA = RandomUniformDistribution(minAngle, numObservations);

    %% Wedge length
    wedgeLength = [0.1, 10]; % Chosen from sensitivity examples
    wedgeLength = [0.1 50]; % To allow for larger variation in z values
    wL = RandomLoguniformDistribution(wedgeLength, numObservations);

    %% Radius
    % Normalise anyway so the important question is at what point can we
    % consider to be a plane wave? - 50m upper limit chosen from
    % sensitivity examples. Schissler use d_s = 0.25m (~size of human head)
    % to determine if edge can shadow. Therefore can probably assume radius
    % of 0.1m an adequate lower limit

    radius = [0.1 50];
    radiusOne = RandomLoguniformDistribution(radius, numObservations);
    radiusTwo = RandomLoguniformDistribution(radius, numObservations);

    rS = min(radiusOne, radiusTwo);
    rR = max(radiusOne, radiusTwo);

    %% z values
    % Ensure reciprocity as swapping S and R is the same, reflecting in the
    % midpoint is the same and reflecting across the wedge is the same.
    sourcePart = rS ./ (rS + rR);

    % Ensure apex on physical edge but can be visible direct around wedge
%     deltaZ = [0 50];
%     dZ = RandomTriangularDistribution(deltaZ, false, numObservations);
%     z = [-(wL / 2 + sourcePart .* dZ) wL / 2 - sourcePart .* dZ];
%     zOne = RandomUniformDistribution(z, numObservations);
%     zTwo = zOne + dZ;

    % COnstrain z values to top and bottom of wedge
    z = [0 * const wL] - wL / 2;
    zOne = RandomUniformDistribution(z, numObservations);
    zTwo = RandomUniformDistribution(z, numObservations);
    dZ = abs(zTwo - zOne);

    zA = zOne + sourcePart .* dZ;
    for i = 1:numObservations
        if zOne(i) < 0
            zOne(i) = -zOne(i);
            zTwo(i) = -zTwo(i);
            zA(i) = -zA(i);
        end
    end
    zS = zOne;
    zR = zTwo;
    truezS = zS + wL / 2;
    truezR = zR + wL / 2;
    truezA = zA + wL / 2;

    input = [deg2rad(wI), deg2rad(bA), deg2rad(mA), wL, rS, rR, zS, zR];

    geometry.wedgeIndex = wI;
    geometry.bendingAngle = bA;
    geometry.minAngle = mA;
    geometry.wedgeLength = wL;
    geometry.rS = rS;
    geometry.rR = rR;
    geometry.zS = zS;
    geometry.zR = zR;
    geometry.zA = zA;
    geometry.truezS = truezS;
    geometry.truezR = truezR;
    geometry.truezA = truezA;
    geometry.thetaS = geometry.minAngle;
    geometry.thetaR = geometry.minAngle + geometry.bendingAngle;
    input = input';
end

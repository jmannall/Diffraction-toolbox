function geometry = RandomGeometryWedge(numObservations)

    % Geometry template
    gtemplate.wedgeIndex = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);

    epsilon = 1e-3;

    %% Wedge index
    numTriObservations = 0.5 * numObservations;
    numUniObservations = numObservations - numTriObservations;

    const = ones(numObservations, 1);

    wedgeIndex = [181, 360 - epsilon];
    wI = [RandomTriangularDistribution(wedgeIndex, false, numTriObservations); RandomUniformDistribution(wedgeIndex, numUniObservations)];
    
    wedgeIndex = [wI, (360 - epsilon) * const];
    wI = RandomTriangularDistribution(wedgeIndex, true, numObservations);
    
    %% Bending angle and minimum angle
    epsilon = 1e-6;
    bendingAngle = [(180 + epsilon) * const, wI - epsilon];
    bA = RandomTriangularDistribution(bendingAngle, false, numObservations);
    
    epsilon = 1e-2;
    minAngle = [epsilon * const, (wI - bA) / 2];
    mA = RandomUniformDistribution(minAngle, numObservations);

    %% Wedge length
    %wedgeLength = [0.1, 10]; % Chosen from sensitivity examples
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
    %radiusOne = RandomUniformDistribution(radius, numObservations);
    %radiusTwo = RandomUniformDistribution(radius, numObservations);

    r1 = min(radiusOne, radiusTwo);
    r2 = max(radiusOne, radiusTwo);

    %% z values
    % Ensure reciprocity as swapping S and R is the same, reflecting in the
    % midpoint is the same and reflecting across the wedge is the same.
    r = r2 + r1;
    r1Part = r1 ./ r;

    % Ensure apex on physical edge but can be visible direct around wedge
    maxL = 100;
    dZ = sqrt(maxL ^ 2 - r .^ 2);
    deltaZ = [0 * const dZ];
    dZ = RandomTriangularDistribution(deltaZ, false, numObservations);
    %dZ = RandomUniformDistribution(deltaZ, numObservations);
    epsilon = 1e-4;
    apex = [epsilon * const wL - epsilon];
    zA = RandomUniformDistribution(apex, numObservations);
    zOne = [zA - r1Part .* dZ, zA + r1Part .* dZ];
    zTwo = [zOne(:,1) + dZ, zOne(:,2) - dZ];

    select = randi(2, 1, numObservations);
    idx = select == 1;
    z1(idx) = zOne(idx);
    z2(idx) = zTwo(idx);
    idx = select == 2;
    z1(idx) = zOne(idx);
    z2(idx) = zTwo(idx);

    z1 = z1';
    z2 = z2';

    for i = 1:numObservations
        if z1(i) > wL(i) / 2
            z1(i) = wL(i) - z1(i);
            z2(i) = wL(i) - z2(i);
            zA(i) = wL(i) - zA(i);
        end
    end

    geometry.wedgeIndex = wI;
    geometry.bendingAngle = bA;
    geometry.minAngle = mA;
    geometry.wedgeLength = wL;
    geometry.rS = r1;
    geometry.rR = r2;
    geometry.zS = z1;
    geometry.zR = z2;
    geometry.zA = zA;
    geometry.thetaS = geometry.minAngle;
    geometry.thetaR = geometry.minAngle + geometry.bendingAngle;
end

% r1 is the shortest radii.
% z1 corresponds to r1. The z-axis zero lies on the corner nearest to z1 and extends into the positive along the edge.

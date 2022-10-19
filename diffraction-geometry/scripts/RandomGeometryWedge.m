function geometry = RandomGeometryWedge(numObservations)

    % Geometry template
    gtemplate.wedgeIndex = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);

    epsilon = 1e-5;

    numTriObservations = 0.5 * numObservations;
    numUniObservations = numObservations - numTriObservations;

    const = ones(numObservations, 1);

    wedgeIndex = [181, 360];
    w = [RandomTriangularDistribution(wedgeIndex, false, numTriObservations); RandomUniformDistribution(wedgeIndex, numUniObservations)];
    
    wedgeIndex = [w, 360 * const];
    w = RandomTriangularDistribution(wedgeIndex, true, numObservations);
    
    bendingAngle = [(180 + epsilon) * const, w];
    bA = RandomTriangularDistribution(bendingAngle, false, numObservations);
    
    minAngle = [epsilon * const, (w - bA) / 2];
    mA = RandomUniformDistribution(minAngle, numObservations);

%     minAngle = [zeros(numObservations, 1), w / 2 - 90];
%     mA = RandomTriangularDistribution(minAngle, false, numObservations);
%     
%     bendingAngle = [180.001 * ones(numObservations, 1), w - 2 * mA];
%     bA = RandomTriangularDistribution(bendingAngle, false, numObservations);

    input = [w, bA, mA];

    geometry.wedgeIndex = input(:,1);
    geometry.bendingAngle = input(:,2);
    geometry.minAngle = input(:,3);
    geometry.source = geometry.minAngle;
    geometry.receiver = geometry.minAngle + geometry.bendingAngle;
end

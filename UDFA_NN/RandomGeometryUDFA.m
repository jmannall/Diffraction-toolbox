% Aim is an representative range of frequency responses (not required to be
% evenly distributed over geometry)

function geometry = RandomGeometryUDFA(numObservations)

    % Geometry template
    gtemplate.wedgeIndex = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);

    %% Wedge index

    epsilon = 1e-2;

    const = ones(numObservations, 1);
    wedgeIndex = [180 + epsilon, 360 - epsilon];

    % pdf = (1 / 4) * (wI - wedgeIndex(1)) .^ 2;
    % Calculate b to ensure the area under the curve equals 1
    cdf = (1 / 12) * (wedgeIndex(2) - wedgeIndex(1)) .^ 3;  % Integral of the pdf
    b = 1 / cdf;
    % New cdf that can be rearranged to make wI the subject and cdf
    % replaced by a uniform distribution.
    %cdf = (b / 12) * (wI - wedgeIndex(1)) .^ 3;

    % wI = nthroot((12 / b) * RandomUniformDistribution([0 1], numObservations), 3) + wedgeIndex(1);
    % b = 12 ./ (wedgeIndex(2) - wedgeIndex(1)) .^ 3;
    wI = (wedgeIndex(2) - wedgeIndex(1)) * nthroot(RandomUniformDistribution([0 1], numObservations), 3) + wedgeIndex(1);

    %% Bending angle and minimum angle

    minAngle = [epsilon * const, (wI - (180 + epsilon)) / 2];
    mA = RandomTriangularDistribution(minAngle, false, numObservations);

    bendingAngle = [epsilon * const, wI - 2 * mA];
    bA = RandomUniformDistribution(bendingAngle, numObservations);
    
    %% Wedge length

    wedgeLength = [0.1 50];
    wL = RandomUniformDistribution(wedgeLength, numObservations);
    % Above 10m roughly inf
    %threshold = 20;
    %wL(wL > threshold) = 5 * wL(wL > threshold) - 4 * threshold;
    
    %% Radius
    % Normalise anyway so the important question is at what point can we
    % consider to be a plane wave? - 50m upper limit chosen from
    % sensitivity examples. Schissler use d_s = 0.25m (~size of human head)
    % to determine if edge can shadow. Therefore can probably assume radius
    % of 0.1m an adequate lower limit

    L = [0.1 100];
    pathLength = RandomUniformDistribution(L, numObservations);

    split = RandomUniformDistribution([epsilon 1 - epsilon], numObservations);

    l = pathLength .* split;
    m = pathLength .* (1 - split);

    apex = [-25 25];
    %apex = [epsilon * const wL-epsilon];

    zA = RandomUniformDistribution(apex, numObservations);
    zA = zA + wL / 2;

    angle = [epsilon 90];
    phii = RandomUniformDistribution(angle, numObservations);

    r1 = l .* sind(phii);
    r2 = m .* sind(phii);

    %% z values
    % Don't need to ensure reciprocity
    z1 = zA - l .* cosd(phii);
    z2 = zA + m .* cosd(phii);

    %% Write struct

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
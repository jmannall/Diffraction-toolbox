function geometry = RandomGeometryWedge_Run2(numObservations)

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
    wedgeIndex = [181, 360 - epsilon];

    % pdf = (1 / 4) * (wI - wedgeIndex(1)) .^ 2;
    % Calculate b to ensure the area under the curve equals 1
    cdf = (1 / 12) * (wedgeIndex(2) - wedgeIndex(1)) .^ 3;  % Integral of the pdf
    b = 1 / cdf;
    % New cdf that can be rearranged to make wI the subject and cdf
    % replaced by a uniform distribution.
    %cdf = (b / 12) * (wI - wedgeIndex(1)) .^ 3;

    wI = nthroot((12 / b) * RandomUniformDistribution([0 1], numObservations), 3) + wedgeIndex(1);

    wI = 180 * (nthroot(RandomUniformDistribution([0 1], numObservations), 3) + 1);
    wI = (wedgeIndex(2) - 180) * nthroot(RandomUniformDistribution([0 1], numObservations), 3) + 180;

    wI(wI < wedgeIndex(1)) = wI(wI < wedgeIndex(1)) + wedgeIndex(1) - 180;
    
%     figure
%     h = histogram(wI, wedgeIndex(1):5:round(wedgeIndex(2)));
%     title('wedgeIndex')
% 
%     n = 2000;
%     w = RandomUniformDistribution(wedgeIndex, n);
%     mA = (w - (180 + epsilon)) / 2;
%     bA = w;
%     area = (1 / 2) * mA .* bA;
% 
%     num = cumsum(round(area * numObservations / sum(area)));
% 
%     num(num > numObservations) = numObservations;
%     num(end) = numObservations;
%     wI = zeros(numObservations, 1);
%     idx = 1:num(1);
%     wI(idx) = w(1);
%     numI(1) = num(1);
%     for i = 2:n
%         idx = num(i - 1):num(i);
%         wI(idx) = w(i);
%         numI(i) = num(i) - num(i - 1);
%     end
% 
%     figure
%     h = histogram(wI, wedgeIndex(1):5:round(wedgeIndex(2)));
%     title('wedgeIndex')
% 
%     figure
%     stem(w, numI)

    %% Bending angle and minimum angle

    minAngle = [epsilon * const, (wI - (180 + epsilon)) / 2];
    mA = RandomTriangularDistribution(minAngle, false, numObservations);

    bendingAngle = [(180 + epsilon) * const, wI - 2 * mA];
    bA = RandomUniformDistribution(bendingAngle, numObservations);
    
    %% Wedge length
    %wedgeLength = [0.1, 10]; % Chosen from sensitivity examples
    wedgeLength = [0.1 50]; % To allow for larger variation in z values
    %wL = RandomLoguniformDistribution(wedgeLength, numObservations);
    wL = RandomUniformDistribution(wedgeLength, numObservations);
    
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

    apex = [epsilon * const wL - epsilon];
    zA = RandomUniformDistribution(apex, numObservations);

    angle = [1 90];
    %phii1 = RandomUniformDistribution(angle, numObservations / 2);
    %phii2 = RandomTriangularDistribution(angle, true, numObservations / 2);
    %phii = [phii1; phii2];
    phii = RandomUniformDistribution(angle, numObservations);

    r1 = l .* sind(phii);
    r2 = m .* sind(phii);

    select = randi(2, numObservations, 1);
    idx = select == 1;
    z1(idx, 1) = zA(idx) - l(idx) .* cosd(phii(idx));
    z2(idx, 1) = zA(idx) + m(idx) .* cosd(phii(idx));
    idx = select == 2;
    z1(idx, 1) = zA(idx) + l(idx) .* cosd(phii(idx));
    z2(idx, 1) = zA(idx) - m(idx) .* cosd(phii(idx));

    %% z values
    % Ensure reciprocity as swapping S and R is the same, reflecting in the
    % midpoint is the same and reflecting across the wedge is the same.
    
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

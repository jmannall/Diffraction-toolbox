function [geometry, ginput] = NNGeometryAugment(gParameters, numObservations)

    constant = ones(numObservations,1);
    epsilon = gParameters.epsilon;

    %% Create wedges

    % Length
    len = length(gParameters.wedgeLength);
    if len == 1
        wedgeLength = gParameters.wedgeLength * constant;
    elseif len == 2
        rWedgeLength = gParameters.wedgeLength;
        wedgeLength = Random(rWedgeLength, numObservations);
        % Transform wedgeLength to weight lower values
        % deltaWedgeLength = rWedgeLength(2) - rWedgeLength(1);
        % wedgeLength = -deltaWedgeLength * sqrt(1 - ((wedgeLength - rangeWedgeLength(1)) / deltaWedgeLength).^2) + rangeWedgeLength(2);
    else
        error('WedgeLength must be a single value or range');
    end

    % Index
    len = length(gParameters.wedgeIndex);
    if len == 1
        wedgeIndex = gParameters.wedgeIndex * constant;
    elseif len == 2
        pd = makedist('Triangular','A',gParameters.wedgeIndex(1).^3,'B',gParameters.wedgeIndex(2).^3,'C',gParameters.wedgeIndex(2).^3);
        wedgeIndex = nthroot(random(pd, numObservations, 1),3);
        %rWedgeIndex = gParameters.wedgeIndex;
        %wedgeIndex = Random(rWedgeIndex, numObservationsgParameters.wedgeIndex(1));
    else
        error('WedgeIndex must be a single value or range');
    end

    %% Create sources and receivers

    len = length(gParameters.thetaS);
    if len == 1
        minAngle = gParameters.thetaS * constant;
    elseif len == 2
        rMinAngle = [gParameters.thetaS(1) * constant, min(wedgeIndex / 2 - epsilon, gParameters.thetaS(2))];
        minAngle = Random(rMinAngle, numObservations);
    else
        error('MinAngle must be a single value or range');
    end

    len = length(gParameters.thetaR);
    if len == 1
        bendingAngle = gParameters.thetaR * constant;
    elseif len == 2
        rBendingAngle = [gParameters.thetaR(1) * constant, min(wedgeIndex - minAngle - epsilon, gParameters.thetaS(2))];
        bendingAngle = Random(rBendingAngle, numObservations);
    else
        error('MinAngle must be a single value or range');
    end


    idx = randi(8,numObservations,1);
    r = gParameters.wedgeIndex;
    store = Random(r, numObservations);
    wedgeIndex(idx == 1) = store(idx == 1);
    wedgeIndex(idx == 2) = store(idx == 2);

    % Shadow boundary
    r = [gParameters.thetaS(1) * constant, (wedgeIndex - 180) / 2 - epsilon];
    store = Random(r, numObservations);
    minAngle(idx == 1) = store(idx == 1);
    r = [-min(5, minAngle), min(5, max(0, wedgeIndex - minAngle - 180))];
    store = Random(r, numObservations);
    bendingAngle(idx == 1) = 180 + store(idx == 1);

    % Specular boundary
    r = [gParameters.thetaS(1), 90 - epsilon];
    store = Random(r, numObservations);
    minAngle(idx == 2) = store(idx == 2);
    r = [-5 * constant, min(5, minAngle)];
    store = Random(r, numObservations);
    bendingAngle(idx == 2) = max(180 - 2 * minAngle(idx == 2) + store(idx == 2), minAngle(idx == 2) + epsilon);

    thetaS = minAngle;
    thetaR = minAngle + bendingAngle;
    % Angle
%     len = length(gParameters.thetaS);
%     if len == 1
%         thetaS = gParameters.thetaS * constant;
%     elseif len == 2
%         rThetaS = [gParameters.thetaS(1) * constant, min(wedgeIndex - epsilon, gParameters.thetaS(2))];
%         thetaS = Random(rThetaS, numObservations);
%     else
%         error('ThetaS must be a single value or range');
%     end
% 
%     thetaR = zeros(size(thetaS));
%     len = length(gParameters.thetaR);
%     if len == 1
%         thetaR = gParameters.thetaR * constant;
%     elseif len == 2
%         for i = 1:length(thetaS)
%             sel = randi(3);
%             if sel == 1 && wedgeIndex(i) - thetaS(i) > 180
%                 rThetaR = [180 - min(5, thetaS(i) / 2), 185];
%             elseif sel == 2 && (wedgeIndex(i) - thetaS(i) > 180 || thetaS(i) > 180)
%                 thetaS(i) = min(wedgeIndex(i) - thetaS(i), thetaS(i));
%                 rThetaR = [thetaS(i) + 175, thetaS(i) + 180 + min(5, thetaS(i) / 2)];
%             else
%                 rThetaR = [gParameters.thetaR(1), min(wedgeIndex(i) - epsilon, gParameters.thetaR(2))];
%             end
%             thetaR(i) = Random(rThetaR, 1);
%         end
%     else
%         error('ThetaR must be a single value or range');
%     end

    % Radius
    len = length(gParameters.radiusS);
    if len == 1
        radiusS = gParameters.radiusS * constant;
    elseif len == 2
        rRadiusS = gParameters.radiusS;
        radiusS = Random(rRadiusS, numObservations);
    else
        error('RadiusS must be a single value or range');
    end

    len = length(gParameters.radiusR);
    if len == 1
        radiusR = gParameters.radiusR * constant;
    elseif len == 2
        rRadiusR = gParameters.radiusS;
        radiusR = Random(rRadiusR, numObservations);
    else
        error('RadiusR must be a single value or range');
    end
    
    % z coordinate
    len = length(gParameters.zS);
    if len == 1
        zS = gParameters.zS * constant;
    elseif len == 2
        rzS = [gParameters.zS(1) * constant, min(wedgeLength - gParameters.zS, gParameters.zS(2))];
        zS = Random(rzS, numObservations);
    else
        error('zS must be a single value or range');
    end

    len = length(gParameters.zR);
    if len == 1
        zR = gParameters.zR * constant;
    elseif len == 2
        rzR = [gParameters.zR(1) * constant, min(wedgeLength - gParameters.zR, gParameters.zR(2))];
        zR = Random(rzR, numObservations);
    else
        error('zR must be a single value or range');
    end

    %% Create NNinput data

    % Checks
    check = wedgeIndex == 180;
    if sum(check) > 0
        wedgeIndex(check) = wedgeIndex(check) + epsilon;
    end
    check = wedgeIndex == 360;
    if sum(check) > 0
        wedgeIndex(check) = wedgeIndex(check) - epsilon;
        thetaR = min(thetaR, wedgeIndex - epsilon);
        thetaS = min(thetaS, wedgeIndex - epsilon);
    end
    check = thetaS == 180;
    if sum(check) > 0
        thetaS(check) = thetaS(check) + epsilon;
    end
    check = wedgeIndex - thetaS == 180;
    if sum(check) > 0
        thetaS(check) = thetaS(check) + epsilon;
    end
    check = thetaR == 180;
    if sum(check) > 0
        thetaR(check) = thetaR(check) + epsilon;
    end
    check = wedgeIndex - thetaR == 180;
    if sum(check) > 0
        thetaR(check) = thetaR(check) + epsilon;
    end
    check = zS == wedgeLength;
    if sum(check) > 0
        zS(check) = zS(check) - epsilon;
    end
    check = zR == wedgeLength;
    if sum(check) > 0
        zR(check) = zR(check) - epsilon;
    end
    check = zS == 0;
    if sum(check) > 0
        zS(check) = epsilon;
    end
    check = zR == 0;
    if sum(check) > 0
        zR(check) = epsilon;
    end

    % Create data
    minAngle = min(thetaS, thetaR);
    bendingAngle = max(thetaS, thetaR) - minAngle;
    visible = double(bendingAngle > 180);
    midPoint = wedgeLength / 2;
    z1 = abs(midPoint - zS);
    z2 = abs(midPoint - zR);
    deltaZ = abs(zR - zS);

    % Bending angle check
    check = bendingAngle == 180;
    if sum(check) > 0
        thetaR(check) = thetaR(check) + epsilon;
        minAngle = min(thetaS, thetaR);
        bendingAngle = max(thetaS, thetaR) - minAngle;
    end

    geometry = struct('wedgeLength', wedgeLength, 'wedgeIndex', wedgeIndex, ...
        'thetaS', thetaS, 'thetaR', thetaR, 'radiusS', radiusS, 'radiusR', radiusR, 'zS', zS, 'zR', zR, 'distribution', 'Augmented');
    
    ginput = struct('wedgeLength', wedgeLength, 'wedgeIndex', wedgeIndex, ...
        'minAngle', minAngle, 'bendingAngle', bendingAngle, 'z1', z1, 'z2', z2, 'deltaZ', deltaZ, 'visible', visible, 'source', minAngle, 'receiver', minAngle + bendingAngle);
end
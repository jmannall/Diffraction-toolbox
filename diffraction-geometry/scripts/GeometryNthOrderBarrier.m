%% Create geometry variable for every combination of the given inputs around an Nth Order barrier.

function geometry = GeometryNthOrderBarrier(barrierRadius, radiusS, radiusR)

    % Geometry template
    gtemplate.barrierRadius = [];
    gtemplate.radiusS = [];
    gtemplate.radiusR = [];

    geometry = repmat(gtemplate, 1, 1);

    [BR, RS, RR] = meshgrid(barrierRadius, radiusS, radiusR);

    bR = reshape(BR,[],1);
    rS = reshape(RS,[],1);
    rR = reshape(RR,[],1);

    store = [bR, rS, rR];

    input = unique(store,'rows');
    
    geometry.barrierRadius = input(:,1);
    geometry.radiusS = input(:,2);
    geometry.radiusR = input(:,3);
end

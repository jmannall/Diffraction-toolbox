%% Create geometry variable for every combination of the given inputs around walls.

function geometry = GeometryWall(wallThickness, bendingAngle, minAngle, radiusR)

    % Geometry template
    gtemplate.wallThickness = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);

    [W, BA, MA] = meshgrid(wallThickness, bendingAngle, minAngle);

    w = reshape(W,[],1);
    bA = reshape(BA,[],1);
    mA = reshape(MA,[],1);

    index =  360 - 2 * asind((wallThickness / 2) / (radiusR + wallThickness / 2)) > mA + bA;
    store = [w, bA, mA];

    input = unique(store(index, :),'rows');
    
    geometry.wallThickness = input(:,1);
    geometry.bendingAngle = input(:,2);
    geometry.minAngle = input(:,3);
    geometry.source = geometry.minAngle;
    geometry.receiver = geometry.minAngle + geometry.bendingAngle;
end

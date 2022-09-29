function geometry = Geometry(wedgeIndex, bendingAngle, minAngle)

    % Geometry template
    gtemplate.wedge = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);

    [W, BA, MA] = meshgrid(wedgeIndex, bendingAngle, minAngle);

    w = reshape(W,[],1);
    bA = reshape(BA,[],1);
    mA = reshape(MA,[],1);

    source = mA;
    receiver = mA + bA;

    index = w - receiver >= mA;
    store = [w, bA, mA];

    input = unique(store(index, :),'rows');
    
    geometry.wedgeIndex = input(:,1);
    geometry.bendingAngle = input(:,2);
    geometry.minAngle = input(:,3);
    geometry.source = geometry.minAngle;
    geometry.receiver = geometry.minAngle + geometry.bendingAngle;
end

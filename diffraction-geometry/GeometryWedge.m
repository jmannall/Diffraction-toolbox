%% Create geometry variable for every combination of the given inputs around wedges.

function geometry = GeometryWedge(wedgeIndex, bendingAngle, minAngle, reciprocity)

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

    if reciprocity
        index = w - bA >= 2 * mA;
    else
        index = w > mA + bA;
    end
    store = [w, bA, mA];

    input = unique(store(index, :),'rows');
    
    geometry.wedgeIndex = input(:,1);
    geometry.bendingAngle = input(:,2);
    geometry.minAngle = input(:,3);
    geometry.source = geometry.minAngle;
    geometry.receiver = geometry.minAngle + geometry.bendingAngle;
end

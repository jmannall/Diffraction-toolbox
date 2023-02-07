%% Create geometry variable for every combination of the given inputs around wedges.

function geometry = GeometryWedge(wedgeIndex, bendingAngle, minAngle, reciprocity, shadowOnly)

    % Geometry template
    gtemplate.wedgeIndex = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);

    [W, BA, MA] = meshgrid(wedgeIndex, bendingAngle, minAngle);

    w = reshape(W,[],1);
    bA = reshape(BA,[],1);
    mA = reshape(MA,[],1);

    index = mA + bA == 180;
    bA(index) = bA(index) + 0.001;

    if reciprocity
        index = w - bA >= 2 * mA;
    else
        index = w > mA + bA;
    end

    if shadowOnly
        indexTemp = bA > 180;
        index = index & indexTemp;
    end

    store = [w, bA, mA];

    input = unique(store(index, :),'rows');
    

    geometry.wedgeIndex = input(:,1);
    geometry.bendingAngle = input(:,2);
    geometry.minAngle = min(input(:,3), input(:,1) - (input(:,3) + input(:,2)));
    geometry.source = input(:,3);
    geometry.receiver = input(:,3) + input(:,2);
end

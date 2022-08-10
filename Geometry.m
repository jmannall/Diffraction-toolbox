
function geometry = Geometry(step, shadowZone)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create geometry function

    % Input data
    wedgeIndex = 180:step:360;
    numWedges = length(wedgeIndex);
    
    % Geometry template
    gtemplate.wedge = [];
    gtemplate.source = [];
    gtemplate.receiver = [];
    gtemplate.bendingAngle = [];
    gtemplate.minAngle = [];
    
    geometry = repmat(gtemplate, 1, 1);
    count = 0;
        
    if shadowZone
        % Shadow zone only
        for i = 1:numWedges
            numSources = (wedgeIndex(i) - 180 + step) / step;
            source = 0:step:(wedgeIndex(i) - 180);
            for j = 1:numSources
                numReceivers = (wedgeIndex(i) - (source(j) + 180) + step) / step;
                receiver = (source(j) + 180):step:wedgeIndex(i);
                for k = 1:numReceivers
                    count = count + 1;
                    geometry(count).wedge = wedgeIndex(i);
                    geometry(count).source = source(j);
                    geometry(count).receiver = receiver(k);
                    geometry(count).bendingAngle = abs(receiver(k) - source(j));
                    geometry(count).minAngle = min(source(j), wedgeIndex(i) - receiver(k));
                end
            end
        end
    else
        % All cases
        for i = 1:numWedges
            numSources = wedgeIndex(i) / step;
            source = 0:step:(wedgeIndex(i) - step);
            for j = 1:numSources
                numReceivers = (wedgeIndex(i) - source(j)) / step;
                receiver = (source(j) + step):step:wedgeIndex(i);
                for k = 1:numReceivers
                    count = count + 1;
                    geometry(count).wedge = wedgeIndex(i);
                    geometry(count).source = source(j);
                    geometry(count).receiver = receiver(k);
                    geometry(count).bendingAngle = abs(receiver(k) - source(j));
                    geometry(count).minAngle = min(source(j), wedgeIndex(i) - receiver(k));
                end
            end
        end
    end
end


function geometry = Geometry(step, shadowZone, minw, maxw)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create geometry function

    check = rem((maxw - minw), step);
    check = 0;

    if check > 0
        disp('Step must be an multiple of the wedge range');
        geometry = 0;
        return
    end

    % Input data
    wedgeIndex = minw:step:maxw;
    numWedges = length(wedgeIndex);
    minsr = step / 2;
    %minsr = 0;
    
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
            %numSources = (wedgeIndex(i) - 180 + minsr) / step;
            source = 0:step:(wedgeIndex(i) / 2);
            numSources = length(source);
            for j = 1:numSources
                %numReceivers = (wedgeIndex(i) - (source(j) + 180) + minsr) / step;
                receiver = (source(j) + 180):step:wedgeIndex(i);
                numReceivers = length(receiver);
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
            %numSources = (wedgeIndex(i) - minsr) / step;
            source = 0:step:(wedgeIndex(i) / 2);
            numSources = length(source);
            for j = 1:numSources
                numReceivers = (wedgeIndex(i) - source(j) - minsr) / step;
                receiver = (source(j) + step):step:wedgeIndex(i);
                numReceivers = length(receiver);
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

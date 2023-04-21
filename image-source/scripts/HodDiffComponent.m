function [hodDiff, plot] = HodDiffComponent(room, plot, diffOrder)
    
    edgeStore = 1:room.numEdges;
    edges = edgeStore;
    for i = 1:diffOrder
        edges = [edges, edgeStore];
    end
    
    hodDiff = cell(diffOrder, 1);
    for j = 2:diffOrder
        diff.valid = false;
        idx = 0;

        tic
        disp(['Diffraction order: ', num2str(j)])
        paths = unique(nchoosek(edges,j), "rows");
        numPaths = size(paths, 1);
        for i = 1:numPaths
            path = paths(i,:);
            if room.sourceCanSeeEdge(path(1)) && room.receiverCanSeeEdge(path(j))
                valid = true;
                n = 1;
                while valid && n < j
                    valid = room.edgeCanSeeEdge(path(n), path(n + 1));
                    n = n + 1;
                end
                if valid
                    diff.valid = true;
                    idx = idx + 1;

                    diff.edges(idx,:) = path;
                    plot.hodDiff{j}{idx} = [room.source; room.edgeMidPoints(path,:); room.receiver];
                end
            end
        end
        toc
        hodDiff{j} = diff;
        clear diff
    end
end
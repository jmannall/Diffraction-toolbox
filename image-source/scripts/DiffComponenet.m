function [diff, plot] = DiffComponenet(room, plot)

    disp('Diffraction order: 1')
    idx = 0;
    %diff.valid = false(room.numEdges, 1);
    diff.valid = false;

    tic
    for i = 1:room.numEdges
        if room.receiverCanSeeEdge(i) && room.sourceCanSeeEdge(i)
            diff.valid = true;
            idx = idx + 1;

            baseCorner = room.corners(room.edgeCorners(i,1),:);
            topCorner = room.corners(room.edgeCorners(i,2),:);
            diff.zW(idx) = PathLength(baseCorner, topCorner);
            d = NormaliseVector(topCorner - baseCorner);
    
            AP = room.source - baseCorner;
            diff.rS(idx) = norm(cross(AP, d));
            diff.zS(idx) = norm(AP) * dot(NormaliseVector(AP), d);
    
            AP = room.receiver - baseCorner;
            diff.rR(idx) = norm(cross(AP, d));
            diff.zR(idx,:) = norm(AP) * dot(NormaliseVector(AP), d);
    
            diff.zA(idx,:) = baseCorner + CalculateApex(diff.rS(idx), diff.rR(idx), diff.zS(idx), diff.zR(idx), diff.zW(idx), true);
    
            kS = NormaliseVector(room.source - diff.zA(idx,:));
            diff.dS(idx) = PathLength(diff.zA(idx,:), room.source);
    
            kR = NormaliseVector(room.receiver - diff.zA(idx,:));
            diff.dR(idx) = PathLength(diff.zA(idx,:), room.receiver);
    
            tS = acosd(dot(kS, room.edgeNormals(i,:)));
            rotS = sign(dot(cross(kS, room.edgeNormals(i,:)), d));
    
            rR = acosd(dot(kR, room.edgeNormals(i,:)));
            rotR = sign(dot(cross(kR, room.edgeNormals(i,:)), d));
    
            diff.tW(idx) = room.thetaW(i);
            halfTW = diff.tW(idx) / 2;
            if rotS == rotR
                diff.bA(idx) = abs(tS - rR);
                if tS < rR
                    diff.tS(idx) = halfTW + tS;
                    diff.tR(idx) = halfTW + rR;
                else
                    diff.tS(idx) = halfTW - tS;
                    diff.tR(idx) = halfTW - rR;
                end
            else
                diff.bA(idx) = tS + rR;
                diff.tS(idx) = halfTW - tS;
                diff.tR(idx) = halfTW + rR;
            end

            diff.edges(idx) = i;
            plot.diff{idx} = [room.source; diff.zA(idx,:); room.receiver];
        end
    end
    toc
end
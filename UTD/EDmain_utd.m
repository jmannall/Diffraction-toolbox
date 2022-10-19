%% UTD model

function [dirmag, diffmag, directlen, pathlen] = EDmain_utd(verts, tris, corners, ps, pk, edgedata, firstorderpathdata, freqs, data, always, xtd_oob, filehandlingparameters)
    [EDversionnumber,lastsavedate,lastsavetime] = EDgetversion;

    Inputdatastruct = struct('corners', verts, 'tris', tris, 'S', ps, 'R', pk, 'data', data, 'EDversionnumber',EDversionnumber);
    EDinputdatahash = DataHash(Inputdatastruct);
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_utd',EDinputdatahash);

    if foundmatch == 1
        eval(['load ''',existingfilename,'''']);
        disp(['   Recycled ',existingfilename]);
    else
        activeedges = find(firstorderpathdata.edgeisactive);
        numactiveedges = length(activeedges);
        numreceivers = size(pk, 1);
        numsources = size(ps, 1);
        [diffmag, dirmag, dir, xm_diff] = deal(zeros(length(freqs), numreceivers));
        [diffir, dirir] = deal(zeros(2 * length(freqs), numreceivers));
        directlen = 1e4 * ones(numreceivers, 1);
        pathlen = 1e4 * ones(numreceivers, numactiveedges);
        tau = 2 * pi;
        m = 1;
        for n = 1:numreceivers
            if numsources ~= 1
                m = n;
            end
            for i = 1:numactiveedges
                edgeverts = corners(edgedata.edgecorners(activeedges(i,:),:),:);
                [ri, rd, zi, zd, li, ld, phii, apex, inedge, zedgelo, zedgehi, thetai, thetad, v] = EDedgeparams(edgeverts, ps(m,:), pk(n,:), edgedata, activeedges(i,:));
                if (thetad - thetai < pi && ~always)
                    fprintf('Not in shadow region');
                end
                if (~inedge && ~xtd_oob)
                    fprintf('Diffraction apex point out of bounds of edge')
                end
                if firstorderpathdata.diffpaths(n,m,activeedges(i)) == 1 && inedge
                    for x = 1:length(freqs)
                        xm = EDutd_core(thetai, li, thetad, ld, v, phii, tau * freqs(x) / data.c, false);
                        if (~always)
                            xm = xm * 2;
                        end
                        xm_diff(x,n) = xm;
                    end
                    pathlen(n,i) = li + ld;
                    diffmag(:,n) = diffmag(:,n) + xm_diff(:,n);
                elseif firstorderpathdata.diffpaths(n,m,activeedges(i)) == 1 && ~inedge
                    pathlen(n,i) = 1000;
                end
            end
            if (~EDsegintersectanytriangle(verts, tris, ps(m,:), pk(n,:)))
                %fprintf('direct path open')
                directlen(n,1) = norm(pk(n,:) - ps(m,:));
                    for x = 1:length(freqs)
                        dir(x,n) = exp(-1i * tau * freqs(x) * directlen(n,1) / data.c) / directlen(n,1);
                    end
                    dirmag(:,n) = dir(:,n);
            end
        end
        desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_utd.mat'];
        eval(['save ''',desiredname,''' dirmag diffmag directlen pathlen EDinputdatahash']);
    end
end
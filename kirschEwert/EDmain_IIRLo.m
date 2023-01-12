function [irtot, pathlenstore] = EDmain_IIRLo(corners, ps, pk, edgedata, firstorderpathdata, data, xtd_oob, filehandlingparameters)
    [EDversionnumber,lastsavedate,lastsavetime] = EDgetversion;

    Inputdatastruct = struct('corners', corners, 'firstorderpathdata', firstorderpathdata, 'S', ps, 'R', pk, 'data', data, 'EDversionnumber',EDversionnumber);
    EDinputdatahash = DataHash(Inputdatastruct);
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_iirlo',EDinputdatahash);

    if foundmatch == 1
        eval(['load ''',existingfilename,'''']);
        disp(['   Recycled ',existingfilename]);
    else
        activeedges = find(firstorderpathdata.edgeisactive);
        numactiveedges = length(activeedges);
        numrecievers = size(pk, 1);
        numsources = size(ps, 1);
        irtot = zeros(data.irlen, numactiveedges, numrecievers);
        filterCoeff = struct('a0', 0, 'a1', 0, 'b0', 0, 'b1', 0);
        T = 1 / data.fs;
        pathlenstore = zeros(numactiveedges, numrecievers);
        m = 1;
        for n = 1:numrecievers
            direct = false;
            count = true;
            if numsources ~= 1
                m = n;
            end
            for i = 1:numactiveedges
                edgeverts = corners(edgedata.edgecorners(activeedges(i,:),:),:);
                [ri, rd, zi, zd, li, ld, phii, apex, inedge, zedgelo, zedgehi, thetai, thetad, v] = EDedgeparams(edgeverts, ps(m,:), pk(n,:), edgedata, activeedges(i,:));
                if (~inedge && ~xtd_oob)
                    fprintf('Diffraction apex point out of bounds of edge')
                end
                if (thetad - thetai < pi)
                    [d, pathlen] = deal(norm(pk(n,:) - ps(m,:)));
                    direct = true;
                else
                    d = (2 * li * ld) / (li + ld);
                    pathlen = li + ld;
                end
                if firstorderpathdata.diffpaths(n,m,activeedges(i)) == 1 && inedge && ((inedge || direct) && count)
                    thetaw = pi / v;
                    thetab = max(0.0001, thetad - (pi + thetai));
                    thetamin = min(thetai, max(0, thetaw - thetab - thetai - pi));
                    thetabmax = thetaw - thetamin - pi;
                    mp = 1 - 0.75 * tanh(1 / (2 * thetab)) * sqrt(tanh(2 * thetamin));
                    mw = (1 - 0.75 * (thetab / thetabmax) * sqrt(sin(-thetaw / 2)))^(-1);
                    fc = (data.c * real(mw) * mp) / (3 * pi * d * (1 - cos(thetab)) * sin(phii) ^ 2);
                    %% LPF
                    f0 = 1.11 * fc * (15.6 / fc) ^ 0.141;
                    LPFfilt = EDlpf(f0, T, filterCoeff);
    
                    %% Hsh
                    actual = (log2(20000) - log2(f0)) * 6;
                    target = (log2(20000) - log2(fc)) * 3;
                    
                    fsh =   209 * fc * (15.6 / fc) ^ 0.827;
                    G = max(0, actual - target);
                    Hshfilt = EDhsh(fsh, G, T, filterCoeff);
    
                    %% Attenuation filter
                    num = 1;
                    den = pathlen;
                    ALLfilt = dfilt.df1(num, den);
    
                    %% Combined filter
                    totfilt = dfilt.cascade(LPFfilt, Hshfilt, ALLfilt);
                    x = [1, zeros(1, data.irlen - 1)];
                    irtot(:,i,n) = transpose(filter(totfilt, x));
                    pathlenstore(i,n) = pathlen;
                    if direct
                        count = false;
                    end
                end
            end
        end
        [irtot, pathlenstore] = EDtrimedges(irtot, pathlenstore);
        desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_iirlo.mat'];
        eval(['save ''',desiredname,''' irtot pathlenstore EDinputdatahash']);
    end
end

function GenerateHOD(idx)

    numExamples = 20;
    nStart = numExamples * idx + 1;
    nEnd = numExamples * (idx + 1);

    load('HODdata.mat')
    
    for i = nStart:nEnd
        wI(i,:) = data(i).wedgeIndex;
        thetaS(i) = data(i).thetaS;
        thetaR(i) = data(i).thetaR;
        bA(i,:) = deg2rad([wI(i,1) - thetaS(i), wI(i,2:numEdges - 1), thetaR(i)]);
    
        [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid, vReceiver] = CreateNthOrderPathData(wI(i,:), thetaS(i,:), thetaR(i,:), rS(i,:), rR(i,:), W(i,:), h);
    
        source(:,3) = zS(i);
        receiver(:,3) = zR(i);
        controlparameters.fs = 2 * fs;
        controlparameters.nfft = 2 * nfft;
        [~, tfmagStore, ~, fvecBtm, ~] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot);
        tfmag.Btm = tfmagStore.diff2;
        tfmagN.Btm(:,i) = CreateNBandMagnitude(tfmag.Btm, fidx);
    
        if zS(i) < zR(i)
            z = zS(i);
        else
            z = zR(i);
        end
        dZ = zR(i) - zS(i);
    
        apex(1,3) = zS(i) + (rS(i)) * dZ / (rS(i) + W(i) + rR(i));
        apex(2,3) = zS(i) + (rS(i) + W(i)) * dZ / (rS(i) + W(i) + rR(i));
    
        if max(apex(:,3)) > h || min(apex(:,3)) < 0
            disp('Apex outside of edge')
        end
    
        % Create virtual sources and receivers
        vSource = [source; apex(1:numEdges - 1,:)];
        vReceiver = [apex(2:numEdges,:); receiver];
    
        rSE = [rS(i), W(i)];
        rRE = [rR(i), fliplr(W(i))];
    
        cumRs = cumsum(rSE)';
        cumRr = fliplr(cumsum(rRE))';
    
        vCorners = [corners(2:numEdges + 1,1:2), vSource(:,3)];
        vector = vSource - vCorners;
        vector(:,3) = 0;
        vector = cumRs .* vector ./ vecnorm(vector, 2, 2);
        vSource = vCorners + vector;
        vSource(:,3) = zS(i);
    
        vector = vReceiver - vCorners;
        vector(:,3) = 0;
        vector = cumRr .* vector ./ vecnorm(vector, 2, 2);
        vReceiver = vCorners + vector;
        vReceiver(:,3) = zR(i);
    
        dZ = cumsum(abs([zS(i) - apex(1,3); apex(1,3) - apex(2,3)]));
        dS = sqrt(cumRs .^ 2 + dZ .^ 2);
        dZ = flipud(cumsum(abs([apex(2,3) - zR(i); apex(1,3) - apex(2,3)])));
        dR = sqrt(cumRr .^ 2 + dZ .^ 2);
        [tfmag.BtmE, ~, tf.BtmE] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, dS, dR, L(i));
        
        % Create virtual sources and receivers
        vSource = [source; apex(1:numEdges - 1,:)];
        vReceiver = [apex(2:numEdges,:); receiver];
    
        dSA = [sqrt(rS(i) ^ 2 + (zS(i) - apex(1,3)) ^ 2), sqrt(W(i) ^ 2 + (apex(1,3) - apex(2,3)))];
        dRA = [sqrt(W(i) ^ 2 + (apex(1,3) - apex(2,3)) ^ 2), sqrt(rR(i) ^ 2 + (zR(i) - apex(2,3)))];
    
        [tfmag.BtmA, ~, tf.BtmA] = ProcessBTMDaisyChain(vSource, vReceiver, corners, planeCorners, controlparameters, numEdges, dSA, dRA, L(i));
        
        %% Interpolated
    
        [tfmag.BtmIE, tf.BtmIE] = deal(zeros(nfft, numEdges + 1));
    
        disp('BTM Ext')
    
        [tfmag.BtmIE(:,1), ~, tf.BtmIE(:,1)] = SingleWedgeInterpolated(h, wI(i,1), epsilon, wI(i,1) - thetaS(i), rS(i), W(i) + rR(i), zS(i), zR(i), controlparameters, false);
        [tfmag.BtmIE(:,2), ~, tf.BtmIE(:,2)] = SingleWedgeInterpolated(h, wI(i,2), epsilon, thetaR(i), rS(i) + W(i), rR(i), zS(i), zR(i), controlparameters, false);
    
    end
end
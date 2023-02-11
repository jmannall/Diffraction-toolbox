function [tfmagOut, fvec, tfcomplex, ir] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot)
    
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    c = controlparameters.c;
    L = 0;
    if isfield(controlparameters, 'Rstart')
        controlparameters.Rstart = controlparameters.Rstart - 1 / fs * c;
        L = controlparameters.Rstart;
        if L < 10
            controlparameters = rmfield(controlparameters, 'Rstart');
        end
    end

    numReceivers = length(thetaR);
    idx = thetaR - thetaS > 180;
    input = [1; zeros(1e4, 1)];

    %% In Shadow

    if sum(idx) > 0
        if L < 10
            zA = CalculateApex(radiusS, radiusR, zS, zR, wedgeLength, true);
            zA = zA(:,3);
    
            pathLength = sqrt(radiusS .^ 2 + (abs(zA - zS) .^ 2)) + sqrt(radiusR .^ 2 + (abs(zA - zR) .^ 2));
            [~, dirIr] = DelayLine(input, pathLength(idx), 3, true(size(pathLength(idx))), c, fs);
            dirRef = IrToTf(dirIr, nfft);
    
            scale = pathLength(idx)';
        else
            dirRef = zeros(nfft / 2, 1);
            scale = 1;
        end
    
        % Generate true response
        [~, tfmag, ~, fvec, tfcomplexBtm] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR(idx), radiusS, radiusR(idx), zS, zR(idx), controlparameters, createPlot);
    
        % Generate reference boundary responses
        epsilon = 1e-10;
        [~, tfmagDiffRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, (thetaS + 180 + epsilon) * ones(size(radiusR(idx))), radiusS, radiusR(idx), zS, zR(idx), controlparameters, createPlot);
    
        % Create scaled response
        diffRef = tfmagDiffRef.diff1;
        shift = dirRef - diffRef;
        scaledResponse = tfmag.diff1 + shift;
        truth = 0;
    
        % Interpolate between the true and scaled responses
        i = min(1, (thetaR(idx) - thetaS - 180) / (wedgeIndex - thetaS - 180 - truth))';
        tfmagOut(:,idx) = (1 - i) .* scaledResponse + i .* tfmag.diff1;
        
        % Adjust magnitude of tfcomplex
        phase = angle(tfcomplexBtm.diff1);
        tfcomplex(:,idx) = scale .* PolarToComplex(10.^(tfmagOut(:,idx) / 20), phase);
        tfmagOut(:,idx) = mag2db(abs(tfcomplex(:,idx)));
    end

    %% Direct

    if sum(~idx) > 0
        source = [radiusS * cosd(thetaS) radiusS * sind(thetaS) zS];
        receiver = [radiusR .* cosd(thetaR) radiusR .* sind(thetaR) zR];
        
        pathLength = vecnorm(source - receiver, 2, 2);
        [~, ir(:,~idx)] = DelayLine(input, pathLength(~idx), 3, true(size(pathLength(~idx))), c, fs);
        dirIr = ir(:,~idx) .* pathLength(~idx)';
        [tfmagOut(:,~idx), tfcomplex(:,~idx)] = IrToTf(dirIr, nfft);
        fvec = fs/nfft*[0:nfft/2-1];
    else
        ir = 0;
    end


%     if thetaR - thetaS > 180
%         % Removed phase shift for very long paths.
%         if L < 10
%             input = [1; zeros(5, 1)];
%             zA = CalculateApex(radiusS, radiusR, zS, zR, wedgeLength, true);
%             zA = zA(3);
% 
%             pathLength = sqrt(radiusS ^ 2 + (abs(zA - zS) ^ 2)) + sqrt(radiusR ^ 2 + (abs(zA - zR) ^ 2));
%             [~, dirIr] = DelayLine(input, pathLength, 3, 1, c, fs);
%             dirRef = IrToTf(dirIr, nfft);
% 
%             scale = pathLength;
%         else
%             dirRef = zeros(nfft / 2, 1);
%             scale = 1;
%         end
% 
%         % Generate true response
%         [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot);
%         
%         % Generate reference boundary responses
%         [~, tfmagDiffRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 + epsilon, radiusS, radiusR, zS, zR, controlparameters, createPlot);
%         
%         %[~, tfmagDirRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 - epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
%         
%         % Create scaled response
%         diffRef = tfmagDiffRef.diff1;
%         shift = dirRef - diffRef;
%         scaledResponse = tfmag.diff1 + shift;
%         truth = 0;
% 
%         % Interpolate between the true and scaled responses
%         i = min(1, (thetaR - thetaS - 180) / (wedgeIndex - thetaS - 180 - truth));
%         tfmagOut = (1 - i) * scaledResponse + i * tfmag.diff1;
%         
%         % Adjust magnitude of tfcomplex
%         phase = angle(tfcomplex.diff1);
%         tfcomplex = scale * PolarToComplex(10.^(tfmagOut / 20), phase);
%         tfmagOut = mag2db(abs(tfcomplex));
%         ir = 0;
%     else
%         % Remove diffracted sound if not in the shadow zone
%         source = [radiusS * cosd(thetaS) radiusS * sind(thetaS) zS];
%         receiver = [radiusR * cosd(thetaR) radiusR * sind(thetaR) zR];
%         
%         pathLength = vecnorm(source - receiver);
%         input = [1; zeros(11, 1)];
%         [~, ir] = DelayLine(input, pathLength, 12, 1, c, fs);
%         dirIr = ir * pathLength;
%         [tfmagOut, tfcomplex] = IrToTf(dirIr, nfft);
%         fvec = fs/nfft*[0:nfft/2-1];
% %         [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
% %         tfcomplex = tfcomplex.direct;
% %         tfmag = tfmag.direct;
%     end
end

function GenerateAudio(audioFile, sceneIdx)

    disp(audioFile)
    
    height = 1.6;
    % checkPoints = [1.5 7 height
    %     5 7 height
    %     6.2 7 height
    %     6.6 6.6 height
    %     8.1 5.1 height];
    
    checkPoints = [1.5 7 height
        5.5 7 height
        6.4 7 height
        7 6.9 height
        7.6 6.6 height
        7.75 6.2 height
        7.9 5 height];
    
    % checkPoints = [1.5 7 height
    %     5.5 7 height
    %     6.3 6.8 height
    %     6.9 6.6 height
    %     7.1 6.2 height
    %     7.2 5 height];
    checkPointSpeeds = [1.3 1.3 1.3 1.3 1.3 1.3 1.3];
    source = [7.5 1 1.6];
    
    updateRate = 100;
    frameRate = 25;
    cattOffset = [6 6 0];
    
    [receivers, receiverHeading] = RyanPathGeneration(checkPoints, checkPointSpeeds, updateRate, frameRate, source, cattOffset);
    
    source = source - cattOffset;
    receivers = receivers - cattOffset;
    
    %% Input data
    
    fs = 48e3;
    nfft = 4096;
    c = 344;
    controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c, 'difforder', 1, 'saveFiles', 3, 'noDirect', false);
    createPlot = false;
    
    %% Geometry data
    
    wedgeLength = 2.5;
    wedgeIndex = 270;
    numReceivers = length(receivers);
    
    %% Create audio
    
    windowLength = 2 * fs / updateRate;
    audioLength = (numReceivers + 1) / updateRate;
    %audioFile = 'musicDoubleBass';
    [audio, audioFs] = LoopAudio(['sourceAudio' filesep audioFile '.wav'], audioLength);
    audio = resample(audio,fs,audioFs);
    
    %% LR testing
    
    % pathLength = 0.001 * ones(1, numReceivers);
    % validPath = true(size(pathLength));
    % tfmag = [0 -6 -18 -12];
    % 
    % audiowrite(['audio', filesep, 'orginal.wav'], audio, fs);
    % PlotSpectrogramOfWAV(['audio', filesep, 'orginal.wav'], [-70 0], nfft);
    % 
    % %% Delay line LR
    % output = DelayLineLR(audio, pathLength, windowLength, validPath, tfmag, c, fs);
    % 
    % audiowrite(['audio', filesep, 'LRfilter.wav'], output, fs);
    % PlotSpectrogramOfWAV(['audio', filesep, 'LRfilter.wav'], [-70 0], nfft);
    
    %% Calculate radius, theta, z
    
    [rS, thetaS, zS] = CalculateGeometryComponents(source, wedgeIndex);
    [rR, thetaR, zR] = CalculateGeometryComponents(receivers, wedgeIndex);
    
    %% Calculate valid paths
    
    % Direct path
    validPath.dir = thetaR <= 180 + thetaS;
    
    % Diffraction paths
    validPath.edI = thetaR > 180 + thetaS;
    validPath.ed = thetaR > 0;
    
    %% Find apex point
    
    [zA, phii] = CalculateApex(rS, rR, zS, zR, wedgeLength, false);
    
    %% Calculate azimuths and elevations for each path
    
    % Diffraction path
    [azimuth.ed, elevation.ed] = CalculateAzimuthElevation(receiverHeading, receivers, zA);
    
    % Direct path
    [azimuth.dir, elevation.dir] = CalculateAzimuthElevation(receiverHeading, receivers, source);
    
    %% Calculate path lengths
    
    % Direct path
    pathLength.dir = PathLength(receivers, source);
    
    % Diffraction paths
    pathLength.ed = DiffractionPathLength(receivers, source, zA);
    
    %% Create delay lines
    
    [delayedAudio, ir] = DelayLine(audio, pathLength, windowLength, validPath, c, fs);
    
    %% BTM data
    
    [ir.ed, ~, ~, ~, tfBtm.ed] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, createPlot);
    
    ir.ed = ir.ed.diff1;
    
    [idxStart, idxEnd, irLength] = deal(zeros(numReceivers, 1));
    for i = 1:numReceivers
        idxStart(i) = find(ir.ed(:,i), 1, "first");
        idxEnd(i) = find(ir.ed(:,i), 1, "last");
        irLength(i) = idxEnd(i) - idxStart(i) + 1;
    end
    nDir = max(irLength);
    irBtm = zeros(nDir, numReceivers);
    for i = 1:numReceivers
        idx = idxStart(i):idxEnd(i);
        irBtm(1:irLength(i),i) = pathLength.ed(i) * ir.ed(idx,i);
    end
    ir.btm = irBtm;
    
    %% BTM audio
    
    %audioPath.ed = ConvolveIR(audio, ir.ed, windowLength, validPath.ed);
    
    %% UTD audio
    
    tfmag = zeros(numReceivers, 4);
    for i = 1:numReceivers
        [tfmag(i,:), ~, ~] = SingleUTDWedgeInterpolated(thetaS, thetaR(i), rS, rR(i), wedgeIndex, phii(i), controlparameters);
    end
    
    %[audioPath.utd, ir.utd] = DelayLineLR(audio, [pathLength.ed, pathLength.dir], windowLength, validPath.edI, tfmag, c, fs, validPath.NN);
    [~, ir.utd] = DelayLineLR(audio, [pathLength.ed, pathLength.dir], windowLength, validPath.ed, tfmag, c, fs, validPath.edI);
    
    %% NN audio
    
    %loadPath = ['NNSaves', filesep, 'IIR-10000_0001-1-09-099-5.mat'];
    %loadPath = ['NNSaves', filesep, 'iir-10067_0001-1-09-099-5-45.mat'];
    %loadPath = ['NNSaves', filesep, 'iir-6118_0001-1-09-099-5-34.mat'];
    %loadPath = ['NNSaves', filesep, '2_iirW-2096_0001-1-09-099-6-16.mat'];
    loadPath = ['NNSaves', filesep, 'iir-2057_0001-1-09-099-3-25.mat'];
    
    load(loadPath, "net");
    
    outputNN = CreateNNOutput(net, wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR, false);
    
    numFilters = 2;
    [z, p, k] = CreateIIRFromNNOutput(outputNN, numFilters);
    [b, a] = IIRFilterCoefficients(z, p, k, numFilters, numReceivers);
    
    [~, ir.NN] = DelayLineIIRFilter(audio, [pathLength.ed, pathLength.dir], windowLength, validPath.ed, b, a, c, fs, validPath.edI);
    %[audioPath.NN, ir.NN] = DelayLineIIRFilter(audio, [pathLength.NN, pathLength.dir], windowLength, validPath.ed, b.NN, a.NN, c, fs, validPath.NN);
    
    %% Transfer functions
    
    [~, tfcomplex] = IrToTf(ir, nfft);
    
    %% NN Reference
    
    [tfcomplex.NNRef, irNNRef] = CalculateNNRef(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, validPath.dir, pathLength.dir, pathLength.ed, false);
    
    % test = [tfcomplex.NNRef; flipud(conj(tfcomplex.NNRef(2:end,:)))];
    % ir.NNRef = ifft(([tfcomplex.NNRef; flipud(conj(tfcomplex.NNRef(2:end,:)))]), 4096);
    
    %% NN loss
    
    % [targetData, fc, fidx] = CreateFrequencyNBands(mag2db(abs(tfcomplex.NN)), fvec, 12);
    % loss = BiquadLoss(output.NN, targetData, numFilters, nfft, fs, fidx);
    % loss = IIRFilterLoss(output.NN, targetData, numFilters, nfft, fs, fidx);
    
    % idx = 50;
    % figure
    % semilogx(fvec, tfmagTest(:,idx))
    % hold on
    % semilogx(fvec, extractdata(tfmagNN(:,idx)))
    % xlim([20 20e3])
    
    
    %% HRTF
    
    disp('Generate HRTFs')
    HRTFFolder = 'HRTF';
    HRTFFile = 'IRC_1059_R_HRIR';
    load([HRTFFolder filesep HRTFFile])
    
    sourcePosition = [l_hrir_S.azim_v l_hrir_S.elev_v];
    hrtfData = permute(cat(3, l_hrir_S.content_m', r_hrir_S.content_m'), [1,3,2]);
    
    hrtfData = permute(resample(hrtfData, fs, l_hrir_S.sampling_hz), [3,2,1]);
    
    hrtf = CreateHRTF(azimuth, elevation, hrtfData, sourcePosition);
    
    %% BRIR data prep
    
    template = struct('L', [], 'R', []);
    cattTemplate = struct();
    brir = struct('catt', cattTemplate, 'dir', template, 'reverb', cattTemplate, 'ed', template, 'edI', template, 'NN', template);
    
    %% Create BRIR
    
    disp('Create BRIRs')
    for i = 1:numReceivers
        if validPath.ed(i)
            brir.ed.L(:,i) = conv(hrtf.ed.L(:,i), ir.ed(:,i));
            brir.ed.R(:,i) = conv(hrtf.ed.R(:,i), ir.ed(:,i));
        end
        if validPath.edI(i)
            brir.utd.L(:,i) = conv(hrtf.ed.L(:,i), ir.utd(:,i));
            brir.utd.R(:,i) = conv(hrtf.ed.R(:,i), ir.utd(:,i));
            brir.NN.L(:,i) = conv(hrtf.ed.L(:,i), ir.NN(:,i));
            brir.NN.R(:,i) = conv(hrtf.ed.R(:,i), ir.NN(:,i));
        else
            brir.utd.L(:,i) = conv(hrtf.dir.L(:,i), ir.utd(:,i));
            brir.utd.R(:,i) = conv(hrtf.dir.R(:,i), ir.utd(:,i));
            brir.NN.L(:,i) = conv(hrtf.dir.L(:,i), ir.NN(:,i));
            brir.NN.R(:,i) = conv(hrtf.dir.R(:,i), ir.NN(:,i));
        end
    end
    
    shift = ceil(-1.95 / c * fs);   % Account for distance of hrtf measurement
    brir.NN.L = circshift(brir.NN.L, shift, 1);
    brir.NN.R = circshift(brir.NN.L, shift, 1);
    brir.utd.L = circshift(brir.utd.L, shift, 1);
    brir.utd.R = circshift(brir.utd.L, shift, 1);
    
    %% Reverb
    
    cattFolder = 'CATTData';
    CheckFileDir(cattFolder)
    brirSave = 'brir_scenes1to6.mat';
    brirFolder = 'BRIR';
    
    loadPath = [brirFolder, filesep, brirSave];
    output = exist([cd filesep loadPath], "file");
    
    if output == 2
        disp('Load CATT BRIRs from mat')
        load(loadPath, "brir")
    else
        disp('BRIR data is missing')
        return
    end
    
    %% Audio
    
    disp('Create audio')
    audioFilePath = ['audio' filesep];
    
    audioOut.NN = ConvolveStereoIR(audio, brir.NN, windowLength);
    audioOut.utd = ConvolveStereoIR(audio, brir.utd, windowLength);
    
    scene = ['scene_' num2str(sceneIdx)];
    disp(scene)
    saveName = [scene '_' audioFile '_'];
    
    audioOut.reverb.(scene) = ConvolveStereoIR(audio, brir.reverb.(scene), windowLength);
    audiowrite([audioFilePath saveName 'reverb.wav'], [audioOut.reverb.(scene).L audioOut.reverb.(scene).R], fs);
    
    audioOut.catt.(scene) = ConvolveStereoIR(audio, brir.catt.(scene), windowLength);
    audiowrite([audioFilePath saveName 'catt.wav'], [audioOut.catt.(scene).L audioOut.catt.(scene).R], fs);
    
    audiowrite([audioFilePath saveName 'NN.wav'], [audioOut.reverb.(scene).L + audioOut.NN.L audioOut.reverb.(scene).R + audioOut.NN.R], fs);
    audiowrite([audioFilePath saveName 'utd.wav'], [audioOut.reverb.(scene).L + audioOut.utd.L audioOut.reverb.(scene).R + audioOut.utd.R], fs);
    audiowrite([audioFilePath saveName 'anchor.wav'], audioOut.reverb.(scene).L + audioOut.reverb.(scene).R, fs);

end
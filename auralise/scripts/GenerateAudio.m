function GenerateAudio(audioFile, sceneIdx)
    
    %% Create audio input

    disp('Load audio file')
    disp(audioFile)

    fs = 48e3;
    updateRate = 100;
    numReceivers = 593;

    windowLength = 2 * fs / updateRate;
    audioLength = (numReceivers + 1) / updateRate;
    [audio, audioFs] = LoopAudio(['sourceAudio' filesep audioFile '.wav'], audioLength);
    audio = resample(audio,fs,audioFs);

    %% Load BRIRs

    disp('Load BRIRs')
    loadPath = ['BRIR' filesep 'brir_scenes1to6'];
    load(loadPath, 'brir')

    %% Audio
    
    disp('Convolve audio')

    
    audioOut.btm = ConvolveStereoIR(audio, brir.ed, windowLength);
    audioOut.utd = ConvolveStereoIR(audio, brir.utd, windowLength);
    audioOut.NN = ConvolveStereoIR(audio, brir.NN, windowLength);
    audioOut.IIR = ConvolveStereoIR(audio, brir.IIR, windowLength);
    
    scene = ['scene_' num2str(sceneIdx)];
    disp(scene)
    
    audioOut.reverb.(scene) = ConvolveStereoIR(audio, brir.reverb.(scene), windowLength);
    audioOut.catt.(scene) = ConvolveStereoIR(audio, brir.catt.(scene), windowLength);
    
    disp('Write audio files')
    audioFolder = 'audio';
    audioFilePath = [audioFolder filesep];
    saveName = [scene '_' audioFile '_'];
    models = fieldnames(audioOut);
    numModels = length(models) - 2;
    
    audiowrite([audioFilePath saveName 'catt.wav'], [audioOut.catt.(scene).L audioOut.catt.(scene).R], fs);
    audiowrite([audioFilePath saveName 'reverb.wav'], [audioOut.reverb.(scene).L audioOut.reverb.(scene).R], fs);
    for i = 1:numModels
        model = models{i};
        audiowrite([audioFilePath saveName model '.wav'], [audioOut.reverb.(scene).L + audioOut.(model).L audioOut.reverb.(scene).R + audioOut.(model).R], fs);
    end
    audiowrite([audioFilePath saveName 'anchor.wav'], audioOut.reverb.(scene).L + audioOut.reverb.(scene).R, fs);
end
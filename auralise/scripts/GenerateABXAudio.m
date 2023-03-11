function GenerateABXAudio(audioFile, sceneIdx, positionIdx, brir)

    %% Create audio input

    disp('Load audio file')
    disp(audioFile)

    fs = 48e3;

    [audio, audioFs] = audioread(['sourceAudio' filesep audioFile '.wav']);
    audio = resample(audio,fs,audioFs);

    %% Load BRIRs

%     disp('Load BRIRs')
%     if sceneIdx == 7
%         loadPath = ['BRIR' filesep 'brir_scene7'];
%     else
%         loadPath = ['BRIR' filesep 'brir_scenes1to6'];
%     end
%     load(loadPath, 'brir')

    %% Audio
    
    
    scene = ['scene_' num2str(sceneIdx)];
    if sceneIdx == 6
        scene = ['scene_' num2str(7)];
    end
    disp(scene)

    disp('Write audio files')
    audioFolder = 'audioABX';
    CheckFileDir(audioFolder)
    audioFilePath = [audioFolder filesep audioFile];
    CheckFileDir(audioFilePath)
    loadModels = {'ed', 'utd' 'NN', 'IIR'};
    saveModels = {'btm', 'utd' 'NN', 'IIR'};
    numModels = length(loadModels);
    numPositions = length(positionIdx);

    for i = 1:numPositions
        if sceneIdx == 6
            scene = ['scene_' num2str(7)];
        end
        audioOut.catt.(scene).L = conv(audio, brir.catt.(scene).L(:,positionIdx(i)));
        audioOut.catt.(scene).R = conv(audio, brir.catt.(scene).R(:,positionIdx(i)));

        audioOut.reverb.(scene).L = conv(audio, brir.reverb.(scene).L(:,positionIdx(i)));
        audioOut.reverb.(scene).R = conv(audio, brir.reverb.(scene).R(:,positionIdx(i)));

        if sceneIdx == 6
            scene = ['scene_' num2str(sceneIdx)];
        end
        saveName = [scene '-' num2str(i) '_' audioFile '_'];
        audiowrite([audioFilePath filesep saveName 'catt.wav'], [audioOut.catt.(scene).L audioOut.catt.(scene).R], fs);
        audiowrite([audioFilePath filesep saveName 'reverb.wav'], [audioOut.reverb.(scene).L audioOut.reverb.(scene).R], fs);

        reverbLength = length(audioOut.reverb.(scene).L);
        for j = 1:numModels
            loadModel = loadModels{j};
            if any(brir.(loadModel).L(:,positionIdx(i)))
                model = saveModels{j};
                audioOut.(model).L = conv(audio, brir.(loadModel).L(:,positionIdx(i)));
                audioOut.(model).R = conv(audio, brir.(loadModel).R(:,positionIdx(i)));
                
                irLength = length(audioOut.(model).L);
                audioOut.(model).L = [audioOut.(model).L; zeros(reverbLength - irLength, 1)];
                audioOut.(model).R = [audioOut.(model).R; zeros(reverbLength - irLength, 1)];

                audiowrite([audioFilePath filesep saveName model '.wav'], [audioOut.reverb.(scene).L + audioOut.(model).L audioOut.reverb.(scene).R + audioOut.(model).R], fs);
            end
        end
    end
end
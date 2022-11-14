load 'ReferenceHRTF.mat' hrtfData sourcePosition

hrtfData = permute(double(hrtfData),[2,3,1]);

sourcePosition = sourcePosition(:,[1,2]);

n = 100;
durationPerPosition = 0.2;
%durationPerPosition = 1;
desiredAz = [-120;-60;0;60;120;0;-120;120];
desiredAz = linspace(90,90,n)';
desiredEl = [-90;90;45;0;-45;0;45;45];
desiredEl = linspace(0,180,n)';

desiredPosition = [desiredAz desiredEl];

interpolatedIR  = interpolateHRTF(hrtfData,sourcePosition,desiredPosition);

leftIR = squeeze(interpolatedIR(:,1,:));
rightIR = squeeze(interpolatedIR(:,2,:));

desiredFs = 48e3;
audioFilePath = 'audio\imagination.wav';
[audio, fs] = LoopAudio(audioFilePath, 0.8 * n * durationPerPosition);
audio = 0.8*resample(audio,desiredFs,fs);
audioFilePath = 'audio\imagination48k.wav';
audiowrite(audioFilePath,audio,desiredFs);

audioLength = length(audio);

fileReader = dsp.AudioFileReader(audioFilePath);
deviceWriter = audioDeviceWriter('SampleRate',fileReader.SampleRate);

leftFilter = dsp.FIRFilter('NumeratorSource','Input port');
rightFilter = dsp.FIRFilter('NumeratorSource','Input port');

samplesPerPosition = durationPerPosition*fileReader.SampleRate;
samplesPerPosition = samplesPerPosition - rem(samplesPerPosition,fileReader.SamplesPerFrame);

sourcePositionIndex = 1;
samplesRead = 0;
while ~isDone(fileReader)
    audioIn = fileReader();
    samplesRead = samplesRead + fileReader.SamplesPerFrame;
    
    leftChannel = leftFilter(audioIn,leftIR(sourcePositionIndex,:));
    rightChannel = rightFilter(audioIn,rightIR(sourcePositionIndex,:));
    
    deviceWriter([leftChannel,rightChannel]);
    
    samplesRemaining = audioLength - samplesRead;
    if sourcePositionIndex > 98
        x = 2;
    end
    if samplesRemaining < 1000
        x = 1;
    end
    if mod(samplesRead,samplesPerPosition) == 0
        sourcePositionIndex = sourcePositionIndex + 1;
    end
end

release(deviceWriter)
release(fileReader)
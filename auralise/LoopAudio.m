%% Output looped or trimmed audio of a set length

function [output, fs, sampleLength] = LoopAudio(filePath, timeLength)
    
    [audio, fs] = audioread(filePath);

    audioLength = length(audio);
    sampleLength = floor(timeLength * fs);
    
    numLoops = sampleLength / audioLength;
    
    output = zeros(sampleLength, 1);
    if numLoops < 1
        output = audio(1:sampleLength);
    else
        for i = 1:floor(numLoops)
            idx = (i - 1) * audioLength + 1:i * audioLength;
            output(idx) = audio;
        end
        extra = audioLength * rem(numLoops, 1);
        idx = i * audioLength + 1:i * audioLength + extra;
        output(idx) = audio(1:extra);
    end
end
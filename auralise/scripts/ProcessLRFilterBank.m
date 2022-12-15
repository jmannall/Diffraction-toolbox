%% Process a single pass through a Linkwitz-Riley filter

function [inputBuffer, outputBuffer, output] = ProcessLRFilterBank(input, b, a, inputBuffer, outputBuffer, tfmag, numBands)

    freqResponse = 10 .^ (tfmag / 20);

%     alpf2 = a(:,1:2,1);
%     blpf2 = b(:,1:2,1);
% 
%     ahpf2 = a(:,1:2,3);
%     bhpf2 = b(:,1:2,3);
% 
%     alpf3 = a(:,3:4,3);
%     blpf3 = b(:,3:4,3);
% 
%     ahpf3 = a(:,3:4,4);
%     bhpf3 = b(:,3:4,4);
% 
%     alpf1 = a(:,3:4,1);
%     blpf1 = b(:,3:4,1);
% 
%     ahpf1 = a(:,3:4,2);
%     bhpf1 = b(:,3:4,2);

    [midOutput, bandOutputs] = deal(zeros(1, numBands));

    for i = 1:2:3
        [inputBuffer(:,1:2,i), outputBuffer(:,1:2,i), filterOutput] = ProcessIIRFilter(input, b(:,1:2,i), a(:,1:2,i), inputBuffer(:,1:2,i), outputBuffer(:,1:2,i));
        % Phase compensating all pass section
        [inputBuffer(:,3:4,i), outputBuffer(:,3:4,i), output1] = ProcessIIRFilter(filterOutput, b(:,3:4,i), a(:,3:4,i), inputBuffer(:,3:4,i), outputBuffer(:,3:4,i));
        [inputBuffer(:,5:6,i), outputBuffer(:,5:6,i), output2] = ProcessIIRFilter(filterOutput, b(:,5:6,i), a(:,5:6,i), inputBuffer(:,5:6,i), outputBuffer(:,5:6,i));
        midOutput(1,i:i+1) = output1 + output2;
    end

    for i = 1:numBands
        [inputBuffer(:,7:8,i), outputBuffer(:,7:8,i), bandOutputs(i)] = ProcessIIRFilter(midOutput(i), b(:,7:8,i), a(:,7:8,i), inputBuffer(:,7:8,i), outputBuffer(:,7:8,i));
    end

    output = sum(freqResponse .* bandOutputs);

%     bandOutputs = zeros(1, numBands);
%     for i = 1:numBands
%         [inputBuffer(:,:,i), outputBuffer(:,:,i), bandOutputs(i)] = ProcessIIRFilter(input, b(:,:,i), a(:,:,i), inputBuffer(:,:,i), outputBuffer(:,:,i));
%     end
%     
%     output = sum(freqResponse .* bandOutputs);
end
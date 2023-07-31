function [H, tfmag, tfcomplex] = LinkwitzRiley(f, fs, nfft, gain)
    T = 1 / fs;
    fc = [176, 775, 3408];

    LPF(1,:) = Butterworth(f, fc(1)) .^ 2;
    LPF(2,:) = Butterworth(f, fc(2)) .^ 2;
    LPF(3,:) = Butterworth(f, fc(3)) .^ 2;

    HPF(1,:) = ButterworthHP(f, fc(1)) .^ 2;
    HPF(2,:) = ButterworthHP(f, fc(2)) .^ 2;
    HPF(3,:) = ButterworthHP(f, fc(3)) .^ 2;

    APF(1,:) = (LPF(3,:) + HPF(3,:));
    APF(2,:) = (LPF(1,:) + HPF(1,:));

    H(1,:) = LPF(2,:) .* APF(1,:) .* LPF(1,:);
    H(2,:) = LPF(2,:) .* APF(1,:) .* HPF(1,:);
    H(3,:) = HPF(2,:) .* APF(2,:) .* LPF(3,:);
    H(4,:) = HPF(2,:) .* APF(2,:) .* HPF(3,:);

    [b(:,1,1), a(:,1,1)] = EDlpf2(fc(2), T);
    [b(:,3,1), a(:,3,1)] = EDlpf2(fc(3), T);
    [b(:,5,1), a(:,5,1)] = EDhpf2(fc(3), T);
    [b(:,7,1), a(:,7,1)] = EDlpf2(fc(1), T);

    [b(:,1,2), a(:,1,2)] = EDlpf2(fc(2), T);
    [b(:,3,2), a(:,3,2)] = EDlpf2(fc(3), T);
    [b(:,5,2), a(:,5,2)] = EDhpf2(fc(3), T);
    [b(:,7,2), a(:,7,2)] = EDhpf2(fc(1), T);

    [b(:,1,3), a(:,1,3)] = EDhpf2(fc(2), T);
    [b(:,3,3), a(:,3,3)] = EDlpf2(fc(1), T);
    [b(:,5,3), a(:,5,3)] = EDhpf2(fc(1), T);
    [b(:,7,3), a(:,7,3)] = EDlpf2(fc(3), T);

    [b(:,1,4), a(:,1,4)] = EDhpf2(fc(2), T);
    [b(:,3,4), a(:,3,4)] = EDlpf2(fc(1), T);
    [b(:,5,4), a(:,5,4)] = EDhpf2(fc(1), T);
    [b(:,7,4), a(:,7,4)] = EDhpf2(fc(3), T);

    for i = 1:4
        b(:,2 * i,:) = b(:,2 * i - 1,:);
        a(:,2 * i,:) = a(:,2 * i - 1,:);
    end

    [tfmag, ~, tfcomplex] = CalculateFilterResponseLR(b, a, nfft, fs);

    [inputBuffer, outputBuffer] = InitialiseLRBuffers();
    input = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    for i = 1:10
        [inputBuffer, outputBuffer, output(i)] = ProcessLRFilterBank(input(i), b, a, inputBuffer, outputBuffer, gain, 4);
    end
end
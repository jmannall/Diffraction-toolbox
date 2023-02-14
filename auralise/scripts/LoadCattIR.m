function brir = LoadCattIR(brir, loadPath, file, idx)
    
    maxIrLength = length(brir.L);
    load([loadPath filesep file])
    if idx < 11
        eval(['irL = h_A0_0', num2str(idx - 1), '_BIN_L;']);
        eval(['irR = h_A0_0', num2str(idx - 1), '_BIN_R;']);
        eval(['clear h_A0_0', num2str(idx - 1), '_BIN_L;']);
        eval(['clear h_A0_0', num2str(idx - 1), '_BIN_R;']);
    else
        eval(['irL = h_A0_', num2str(idx - 1), '_BIN_L;']);
        eval(['irR = h_A0_', num2str(idx - 1), '_BIN_R;']);
        eval(['clear h_A0_', num2str(idx - 1), '_BIN_L;']);
        eval(['clear h_A0_', num2str(idx - 1), '_BIN_R;']);
    end
    irLength = length(irL);
    maxIrLength = max(irLength, maxIrLength);
    irL = [irL zeros(1, maxIrLength - irLength)];
    irR = [irR zeros(1, maxIrLength - irLength)];
    brir.L = [brir.L irL'];
    brir.R = [brir.R irR'];
end
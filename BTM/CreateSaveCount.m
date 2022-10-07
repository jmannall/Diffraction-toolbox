function saveCount = CreateSaveCount(result, numInputs, numSaves, filesPerSave)    
    
    count = result(end).i;
    if length(result) == numInputs    % When full result loaded
        saveCount = numSaves + 1;
        disp('Load array from save');
    elseif count == numInputs         % When all results complete but not compiled yet
        saveCount = floor(count / filesPerSave);
        disp('Load previous progress from save');
    else                              % When some results completed
        saveCount = floor(count / filesPerSave) + 1;
        disp('Load previous progress from save');
    end
end
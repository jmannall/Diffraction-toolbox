%% MSE function across logarithmic frequencies

function error = Error(tfinput, tffit, fvec)

    % Distribute frequencies evenly across the audible spectrum
    errorFreq = logspace(log10(20), log10(20000), 100);

    % Create frequency variables
    numfreqs = length(errorFreq);
    freqStore = zeros(1,numfreqs);
    
    % Loop through errorFreq and find closest match in fvec
    for i = 1:numfreqs
        a = find(errorFreq(i) < fvec);
        % Store frequency index
        freqStore(i) = a(1);
    end

    % Mean squared error between input and fit tf.
    error = sum((tfinput(freqStore) - tffit(freqStore)) .^ 2) / numfreqs;
end
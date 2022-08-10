function error = Error(finput, ffit, fvec)
% Not what need. Need JND as it varies with frequency - Measurement and Analysis of Just Noticeable
% Difference of Interaural Level Difference Cue paper sugegsts sensitivity
% decreases above 6kHz
% [spl splfvec] = iso226(60);
% invspl = spl(18) - spl;
errorfreq = logspace(1.3, 4.3, 100);
numfreqs = length(errorfreq);
freqs_store = zeros(1,numfreqs);
error = 0;

for i = 1:numfreqs
    a = find(errorfreq(i) < fvec);
    freqs_store(i) = a(1);
end


error = sum(abs(finput(freqs_store) - ffit(freqs_store))) / numfreqs;

end
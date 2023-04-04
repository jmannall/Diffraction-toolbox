close all
% 
% fs = 48e3;
% nfft = 8192;
% c= 344;
% controlparameters = struct('fs', 2 * fs, 'nfft', 2 * nfft, 'difforder', 1, 'c', c, 'saveFiles', 3, 'noDirect', false);
% 
% [~, tfmagDefault, ~, fvec] = DefaultBTM(controlparameters);
% nBands = 8;
% tfmagDefault = tfmagDefault(1:end / 2,:);
% 
% fvec = fs/nfft*[0:nfft/2-1];
% [~, fc, fidx] = CreateFrequencyNBands(tfmagDefault, fvec, nBands);
% 
% AWeighting = weightingFilter('A-weighting',fs);
% aWeight = freqz(AWeighting, fc, fs);
% aWeight = abs(aWeight)';
% scale = sum(aWeight);
% aWeightNorm = aWeight / scale;
% 
% figure
% semilogx(fc, aWeight)
% hold on
% semilogx(fc, aWeightNorm)
% grid on

alpha = [1 0.5 0.8];

%alpha = alpha / sum(alpha);
Ht = [1 1 1];
Hp = [2 0.5 1.1];

k = length(alpha);

err = abs(20 * log10(Ht ./ Hp));

error = 20 * log10(sum(10 .^ (err / 20) .* alpha) / k)

error = 20 * log10(sum(10 .^ (err / 20) .* alpha))

error = 20 * log10(sum(10 .^ (err / 20)) / k)

error = 20 * log10(sum(10 .^ (err / 20) .* alpha) / sum(alpha))

for i = 1:k
    idx = 2 * i - 1:2 * i;
    x(idx) = i:i + 1;
    HtPlot(idx) = 20 * log10(Ht(i));
    HpPlot(idx) = 20 * log10(Hp(i));
    alphaPlot(idx) = 20 * log10(alpha(i));
end

figure
plot(x, HtPlot)
hold on
plot(x, HpPlot, '--')
plot(x, alphaPlot, '-.')
grid on
ylim([-9 15])
legend('target', 'prediction', 'weighting')

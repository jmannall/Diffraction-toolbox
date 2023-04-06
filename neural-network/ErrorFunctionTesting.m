close all

alpha = [1 0.5 0.8];

%alpha = alpha / sum(alpha);
Ht = [1 1 1];
Hp = [2 0.6 1.1];

k = length(alpha);

err = abs(20 * log10(Ht ./ Hp));

errorAlphaK = 20 * log10(sum(10 .^ (err / 20) .* alpha) / k)

errorAlpha = 20 * log10(sum(10 .^ (err / 20) .* alpha))

errorK = 20 * log10(sum(10 .^ (err / 20)) / k)

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

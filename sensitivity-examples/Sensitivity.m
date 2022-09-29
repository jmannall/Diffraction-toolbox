
wedgeLength = 10;
wedgeIndex = 270;
thetaS = 10;
thetaR = 260;
radiusS = [2,3, 1];
radiusR = [2,1, 3];
zR = 5;
zS = 5;
fs = 48000;

result = struct();
for i = 1:length(radiusR)
    [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS(i), radiusR(i), zS, zR, fs);
end

%% Process data
tfmag = [result.tfmag];
ir = {result.ir};
fvec = [result(1).fvec];
tvec = {result.tvec};

%% Plot
close all

figure
semilogx(fvec, tfmag)
l = legend(string(round(radiusR,2)),'Location','eastoutside');
title(l, 'radiusR')
title(['fs: ', num2str(fs)])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on

figure
for i = 1:length(ir)
    plot(tvec{i}, ir{i})
    hold on
end
hold off
l = legend(string(round(radiusR,2)),'Location','eastoutside');
title(l, 'radiusR')
title(['fs: ', num2str(fs)])
xlabel('Time (s)')
ylabel('Magnitude (dB)')
grid on

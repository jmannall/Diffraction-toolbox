function SensitivityzRzSsame(fs)
    wedgeLength = 10;
    wedgeIndex = 270;
    thetaS = 10;
    thetaR = 260;
    radiusS = 1;
    radiusR = 1;
    zR = -6:3;
    zR(zR == 0) = 0.001;
    zS = zR + 5;
    zS(zS == 0) = 0.001;
    
    result = struct();
    for i = 1:length(zR)
        [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS(i), zR(i), fs);
    end
    
    %% Process data
    tfmag = [result.tfmag];
    ir = {result.ir};
    fvec = [result(1).fvec];
    tvec = {result.tvec};
    
    %% Plot
    figure
    semilogx(fvec, tfmag)
    l = legend(string(round(zR,2)),'Location','eastoutside');
    title(l, 'zR')
    title(['fs: ', num2str(fs)])
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    xlim([20 20000])
    ylim([-60 0])
    grid on
    
    figure
    for i = 1:length(ir)
        plot(tvec{i}, ir{i})
        hold on
    end
    hold off
    l = legend(string(round(zR,2)),'Location','eastoutside');
    title(l, 'zR')
    title(['fs: ', num2str(fs)])
    xlabel('Time (s)')
    ylabel('Magnitude (dB)')
    grid on
end
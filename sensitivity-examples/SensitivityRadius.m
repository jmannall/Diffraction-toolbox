function SensitivityRadius(fs)
    wedgeLength = 10;
    wedgeIndex = 270;
    thetaS = 10;
    thetaR = 260;
    radiusS = 1;
    radiusR = 0.25:0.25:5;
    zS = 5;
    zR = 5;
    
    result = struct();
    for i = 1:length(radiusR)
        [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR(i), zS, zR, fs);
    end
    
    %% Process data
    tfmag = [result.tfmag];
    ir = {result.ir};
    fvec = [result(1).fvec];
    tvec = {result.tvec};
    
    %% Plot    
    figure
    semilogx(fvec, tfmag)
    l = legend(string(radiusR),'Location','eastoutside');
    title(l, 'radiusR')
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
    l = legend(string(radiusR),'Location','eastoutside');
    title(l, 'radiusR')
    title(['fs: ', num2str(fs)])
    xlabel('Time (s)')
    ylabel('Magnitude (dB)')
    grid on
end
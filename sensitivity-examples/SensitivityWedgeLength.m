function SensitivityWedgeLength(fs)
    wedgeLength = logspace(log10(0.1), 2, 20);
    wedgeIndex = 270;
    thetaS = 10;
    thetaR = 260;
    radiusS = 1;
    radiusR = 1;
    zS = wedgeLength / 2;
    zR = wedgeLength / 2;
    
    result = struct();
    for i = 1:length(wedgeLength)
        [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec] = SingleWedge(wedgeLength(i), wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS(i), zR(i), fs);
    end
    
    %% Process data
    tfmag = [result.tfmag];
    ir = {result.ir};
    fvec = [result(1).fvec];
    tvec = {result.tvec};
    
    %% Plot    
    figure
    semilogx(fvec, tfmag)
    l = legend(string(round(wedgeLength,2)),'Location','eastoutside');
    title(l, 'Wedge Length')
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
    l = legend(string(round(wedgeLength,2)),'Location','eastoutside');
    title(l, 'Wedge Length')
    title(['fs: ', num2str(fs)])
    xlabel('Time (s)')
    ylabel('Magnitude (dB)')
    grid on
end
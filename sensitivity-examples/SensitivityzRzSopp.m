function SensitivityzRzSopp(fs)
    wedgeLength = 20;
    wedgeIndex = 270;
    thetaS = 10;
    thetaR = 260;
    radiusS = 1;
    radiusR = 1;
    zR = wedgeLength / 2 + linspace(4, 35, 10);
    zR(zR == 0) = 0.001;
    zS = wedgeLength-zR;
    
    controlparameters = struct('fs', fs, 'nfft', 16384, 'difforder', 1);

    result = struct();
    for i = 1:length(zR)
        [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec, result(i).tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS(i), zR(i), controlparameters, false);
        L = 2 * sqrt(radiusR ^ 2 + (zR(i) - wedgeLength / 2) ^ 2);
        result(i).tfcomplex.diff1 = L .* result(i).tfcomplex.diff1;
    end
    
    %% Process data
    tfcomplex = [result.tfcomplex];
    tfmag = mag2db(abs([tfcomplex.diff1]));
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
        plot(tvec{i}, [ir{i}.diff1])
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
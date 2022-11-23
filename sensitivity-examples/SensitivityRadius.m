function SensitivityRadius(fs)
    wedgeLength = 20;
    wedgeIndex = 270;
    thetaS = 10;
    thetaR = 250;
    radiusS = [10 20 30 40 50];
    radiusR = [10 20 30 40 50];
    zS = 10;
    zR = 10;
    
    controlparameters = struct('fs', fs, 'nfft', 16384, 'difforder', 1, 'saveFiles', 2);

    result = struct();
    for i = 1:length(radiusR)
        controlparameters.Rstart = radiusS(i) + radiusR(i);
        [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec, result(i).tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS(i), radiusR(i), zS, zR, controlparameters, false);
        %result(i).tfcomplex.diff1 = (radiusS(i) + radiusR(i)) .* result(i).tfcomplex.diff1;
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
        plot(tvec{i}, [ir{i}.diff1])
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
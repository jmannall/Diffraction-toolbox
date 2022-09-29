function BTMDataPlot(dVar, iVar1, iVar2, dVarName, gradient, fc)

    toFind = [iVar1, iVar2];
    [~, ~, idx] = unique(toFind, 'rows');
    
    numUnique = max(idx);
    figure
    title('Gradient')
    ylabel('dB/octave')
    xlabel(dVarName)
    xlim([min(dVar) max(dVar)])
    ylim([-5 -2])
    hold on
    for i = 1:numUnique
        index = idx == i;
        plot(dVar(index), gradient(index))
    end    
    hold off

    
    figure
    for i = 1:numUnique
        index = idx == i;
        semilogy(dVar(index), fc(index))
        hold on
    end

    title('Frequency')
    ylabel('Hz')
    ylim([20 2000])
    xlim([min(dVar) max(dVar)])
    xlabel(dVarName)
    hold off

end
%% Plot a spectogram 

function PlotSpectogram(x, f, t, limits, titleText, xLabel)

    if nargin < 6
        xLabel = 'Time (s)';
    end

    % Create figure
    figure
    grid on
    sh = surf(t, f, x);
    c = colorbar;
    view([0 90])
    axis tight
    xlabel(xLabel)
    ylabel('Frequency (Hz)')
    set(gca,'YScale','log')
    set(sh, 'EdgeColor','none')
    clim(limits)
    c.Label.String = 'Magnitude (dB)';
    title(titleText)
end
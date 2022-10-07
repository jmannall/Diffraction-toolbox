%% Plot a spectogram 

function PlotSpectogram(tfcomplex, f, t, limits, titleText, phase, xLabel)

    if nargin < 7
        xLabel = 'Time (s)';
    end
    
    x = mag2db(abs(tfcomplex));
    x = max(x, limits(1));
    position = [50 50 700 600];
    titleText = strrep(titleText, '_', ' ');
    if phase
        y = angle(tfcomplex);
        
        position(3) = 2 * position(3);
        figure("Position", position);
        tl = tiledlayout(1, 2);
        nexttile
        sh = surf(t, f, x);
        grid on
        c = colorbar;
        view([0 90])
        axis tight
        xlabel(xLabel)
        ylabel('Frequency (Hz)')
        xlim([min(t) max(t)])
        set(gca,'YScale','log')
        set(sh, 'EdgeColor','none')
        clim(limits)
        c.Label.String = 'Magnitude (dB)';
        title('Magnitude')

        nexttile
        sh = surf(t, f, y);
        grid on
        c = colorbar;
        view([0 90])
        axis tight
        xlabel(xLabel)
        ylabel('Frequency (Hz)')
        xlim([min(t) max(t)])
        set(gca,'YScale','log')
        set(sh, 'EdgeColor','none')
        clim([-pi pi])
        c.Label.String = 'Phase (rad)';
        title('Phase')

        title(tl, titleText)
    else
        % Create figure
        figure("Position", position);
        sh = surf(t, f, x);
        grid on
        c = colorbar;
        view([0 90])
        axis tight
        xlabel(xLabel)
        ylabel('Frequency (Hz)')
        xlim([min(t) max(t)])
        set(gca,'YScale','log')
        set(sh, 'EdgeColor','none')
        clim(limits)
        c.Label.String = 'Magnitude (dB)';
        title(titleText)
    end

    if ~exist('figures', 'dir')
       mkdir figures
    end
    saveas(gcf, ['figures\', titleText, '.png'])
end
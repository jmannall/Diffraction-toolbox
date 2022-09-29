function plotFilterParameters(geometry, igeometry, variable, fclpf, fchsh, B0, c1Name, c2Name)

    for i = 1:length(geometry)
        index = igeometry == i;
        figure('Units','normalized','Position',[0 0 1 1]);
        sgtitle([c1Name, ': ', num2str(geometry(i,1)), ' ', c2Name, ': ', num2str(geometry(i,2))])
        subplot(1,3,1)
        plot(variable(index),fclpf(index))
        xlim([min(variable) max(variable)])
        ylim([min(fclpf) max(fclpf)]);
        title('fc lpf');

        subplot(1,3,2)
        plot(variable(index),fchsh(index))
        xlim([min(variable) max(variable)])
        ylim([min(fchsh) max(fchsh)]);
        title('fc hsh');

        subplot(1,3,3)
        plot(variable(index),B0(index))
        xlim([min(variable) max(variable)])
        ylim([min(B0) max(B0)]);
        title('B0');

        saveas(gcf,['figures/pso_results/','Parameters_' c1Name, '_', num2str(geometry(i,1)), '_' c2Name, '_', num2str(geometry(i,2)), '.png'])
        close all
    end
end
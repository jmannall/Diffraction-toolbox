function plotFilterPZK(geometry, igeometry, variable, z, p, k, c1Name, c2Name, zPrediction, pPrediction, kPrediction)

    newcolors = [0 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560
        0.4660 0.6740 0.1880];

    for i = 1:length(geometry)
        index = igeometry == i;
        figure
        colororder(newcolors)
        plot(variable(index),z(1,index))
        hold on
        plot(variable(index),z(2,index))
        plot(variable(index),p(1,index))
        plot(variable(index),p(2,index))
        plot(variable(index),k(index))
        % Prediction
        plot(variable(index),zPrediction(1,index),'--')
        plot(variable(index),zPrediction(2,index),'--')
        plot(variable(index),pPrediction(1,index),'--')
        plot(variable(index),pPrediction(2,index),'--')
        plot(variable(index),kPrediction(index),'--')
        hold off
        title([c1Name, ': ', num2str(geometry(i,1)), ' ', c2Name, ': ', num2str(geometry(i,2))])
        legend('z1', 'z2', 'p1', 'p2', 'k')
        xlim([min(variable) max(variable)])
        ylim([-1 1])
    
        saveas(gcf,['figures/pso_results/', c1Name, '_', num2str(geometry(i,1)), '_' c2Name, '_', num2str(geometry(i,2)), '.png'])
        close all
    end
end
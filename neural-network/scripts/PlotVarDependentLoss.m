function PlotVarDependentLoss(data, variable, width, percentiles, titleText, cIdx)

    maxWl = max(variable);
    
    numBins = ceil(maxWl / width);
    [x, meanAv, medianAv] = deal(zeros(1, numBins));
    
    numPercentiles = length(percentiles);
    for i = 1:numBins
        idx = width * (i - 1) < variable & variable <= width * i;
        input = data(idx);
        meanAv(i) = mean(input);
        medianAv(i) = median(input);
        for j = 1:numPercentiles
            p = percentiles(j);
            pTop(j,i) = prctile(input, 100-p);
            pBot(j,i) = prctile(input, p);
        end
        x(i) = (i - 0.5) * width;
    end
    x = [0, x, width * numBins];
    meanAv = [meanAv(1), meanAv, meanAv(end)];
    medianAv = [medianAv(1), medianAv, medianAv(end)];
    pTop = [pTop(:,1), pTop, pTop(:,end)];
    pBot = [pBot(:,1), pBot, pBot(:,end)];
    
    ax = gca;
    alphaValue = 0.15;
    colorValue = ax.ColorOrder(cIdx,:);
    
    for j = 1:numPercentiles
        for i=1:numBins + 1
            Px = [x(i) x(i+1) x(i+1) x(i)];
            Py = [pTop(j,i) pTop(j,i+1) pBot(j,i+1) pBot(j,i)];
            fill(Px,Py,colorValue,'FaceAlpha',alphaValue,'EdgeColor','none');
            hold on
        end
    end

    plot(x, medianAv, 'Color',colorValue)
    plot(x, meanAv, '--', 'Color',colorValue)
    grid on
    ylim([0 6])
    title(titleText)
end
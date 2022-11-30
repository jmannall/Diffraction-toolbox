function x = LoadBayesoptResult(size, filter)
    
    filePath = ['bayesoptResults', filesep, 'BayesoptResult_Size_'];
    load([filePath, num2str(size), '_Filter_', filter], 'lossEst', 'lossObs', 'xEst', 'xObs')
    if lossEst < lossObs
        x = xEst;
    else
        x = xObs;
    end
end
%% Custom ADAM gradient function

function [params, averageGrad, averageSqGrad] = ADAM(params,grad,averageGrad,averageSqGrad,iteration, learnRate, gradDecay, gradSqDecay)   
    epsilon = 1e-8;
    gradDecayT = gradDecay ^ iteration;
    gradSqDecayT = gradSqDecay ^ iteration;

    numCells = numel(params);

    if iteration < 1 || isinteger(iteration)
        disp('Iteration must be an integer greater than zero')
        return
    end
    if iteration == 1
        numCellsX = size(params, 1);
        numCellsY = size(params, 2);
        [averageGrad, averageSqGrad] = deal(cell(numCellsX,numCellsY));
        for i = 1:numCells
            paramsSize = params{i};
            averageGrad{i} = zeros(size(paramsSize));
            averageSqGrad{i} = zeros(size(paramsSize));
        end
    end
    
    for i = 1:numCells
        gradSq = grad{i} .^ 2;
        averageGrad{i} = gradDecay .* averageGrad{i} + (1 - gradDecay) .* grad{i};
        averageSqGrad{i} = gradSqDecay .* averageSqGrad{i} + (1 - gradSqDecay) .* gradSq;
        averageGradVec = averageGrad{i} ./ (1 - gradDecayT);
        averageSqGradVec = averageSqGrad{i} ./ (1 - gradSqDecayT);
        params{i} = params{i} - learnRate .* averageGradVec ./ (sqrt(averageSqGradVec) + epsilon);
    end
end
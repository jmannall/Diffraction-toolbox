%% Train custom neural network

function self = Train(self, Loss)
    start = tic;
    figure
    C = colororder;
    lineLossTrain = animatedline('Color',C(2,:));
    lineLossEpoch = animatedline('Color',C(1,:));
    ylim([0 inf])
    xlabel("Iteration")
    ylabel("Loss")
    grid on
    iteration = 0;
    averageGrad = [];
    averageSqGrad = [];
    gradDecay = 0.9;
    gradSqDecay = 0.999;

    epochLosses = zeros(1, self.numEpochs);
    i = 0;
    epochIncreasing = 0;
    count = 1;
    while i < self.numEpochs
        % Shuffle data.
        i = i + 1;
        idx = randperm(size(self.targetData, 1));
        self.targetData = self.targetData(idx,:);
        self.trainingData = self.trainingData(idx,:);
        for j = 1:self.numIterationsPerEpoch
            idx = (j-1)*self.miniBatchSize+1:j*self.miniBatchSize;
            self.X = self.trainingData(idx,:);
            self.T = self.targetData(idx,:);
            iteration = iteration + 1;
            self = Forward(self);
            % self = UpdateWeights(self, @Delsigmoid);
            [self, averageGrad, averageSqGrad] = UpdateWeightsADAM(self, averageGrad, averageSqGrad, iteration, gradDecay, gradSqDecay, Loss);
            if mod(iteration, 100) == 0
                idx = self.numIterationsPerEpoch * (i - 1) + j;
                loss = self.losses(idx);
                D = duration(0,0,toc(start),'Format',"hh:mm:ss");
                addpoints(lineLossTrain,idx,loss)
                title("Epoch: " + i + ", Elapsed: " + string(D))
                drawnow
                disp(['Loss: ', num2str(loss)]);
            end
        end
        epochLosses(i) = sum(self.losses(((i - 1) * self.numIterationsPerEpoch + 1):i * self.numIterationsPerEpoch)) / self.numIterationsPerEpoch;
        idx = self.numIterationsPerEpoch * i;
        epochLoss = epochLosses(i);
        D = duration(0,0,toc(start),'Format',"hh:mm:ss");
        addpoints(lineLossEpoch,idx,epochLoss)
        title("Epoch: " + i + ", Elapsed: " + string(D))
        drawnow
        disp(['Epoch loss: ', num2str(epochLoss)]);
        avEpochLoss = sum(epochLosses((max(1,i - 9)):i)) / 10;
        if i > 10 && epochLoss > avEpochLoss
            self.learnRate = max(self.learnRate * 0.5, 1e-5);
            count = count + 1;
        end
    end
    if i == self.numEpochs
        disp('Iterations complete')
    else
        disp('Loss no longer decreasing');
    end
end
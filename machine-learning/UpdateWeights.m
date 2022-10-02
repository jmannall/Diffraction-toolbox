%% Simple update weights function for MSE error function

function self = UpdateWeights(self, Delactivation)
    loss = 0.5 * (self.T - self.outputFinal).^2;
    self.losses(end+1) = sum(loss);

    delloss = self.T - self.outputFinal;

    grad01 = self.trainingData' * (((delloss .* Delactivation(self.output)) .* self.weightsOutput') .* Delactivation(self.hidden{1}));
    grad12 = self.hiddenOut{1}' * (delloss .* Delactivation(self.output));
    
    self.weightsLayer{1} = self.weightsLayer{1} + self.learnRate .* grad01;
    self.weightsOutput = self.weightsOutput + self.learnRate .* grad12;

    self.biasLayer{1} = self.biasLayer{1} + self.learnRate .* sum(((delloss .* Delactivation(self.output)) .* self.weightsOutput') .* Delactivation(self.hidden{1}), 1);
    self.biasOutput = self.biasOutput + self.learnRate .* sum(delloss .* Delactivation(self.outputFinal), 1);
end

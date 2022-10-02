function [self, averageGrad, averageSqGrad] = UpdateWeightsADAM_Temp(self, averageGrad, averageSqGrad, iteration, gradDecay, gradSqDecay, Loss)
    loss = 0.5 * (self.outputFinal - self.targetData).^2;
    delloss = self.outputFinal - self.targetData;

%     [loss, delloss] = dlfeval(Loss, dlarray(self.outputFinal), self.targetData);
    self.losses(end+1) = sum(loss, "all");
    
%     nodeValue12 = delloss .* self.Delactivation(self.output);
%     gradw12 = self.hiddenOut{1}' * nodeValue12;
%     gradb12 = sum(nodeValue12, 1);
%     nodeValue01 = nodeValue12 * self.weightsOutput' .* self.Delactivation(self.hidden{1});
%     gradw01 = self.trainingData' * nodeValue01;
%     gradb01 = sum(nodeValue01, 1);
    gradw01 = self.trainingData' * (((delloss .* self.Delactivation(self.output)) * self.weightsOutput') .* self.Delactivation(self.hidden{1}));
    gradw12 = self.hiddenOut{1}' * (delloss .* self.Delactivation(self.output));

    gradb01 = sum(((delloss .* self.Delactivation(self.output)) * self.weightsOutput') .* self.Delactivation(self.hidden{1}), 1);
    gradb12 = sum(delloss .* self.Delactivation(self.outputFinal), 1); %Should
    %this one be self.output?

    params = {self.weightsLayer{1}, self.biasLayer{1}; self.weightsOutput, self.biasOutput};
    grad = {gradw01, gradb01; gradw12, gradb12};

    [params, averageGrad, averageSqGrad] = ADAM(params,grad,averageGrad,averageSqGrad,iteration,self.learnRate,gradDecay,gradSqDecay);
    self.weightsLayer{1} = params{1,1};
    self.weightsOutput = params{2,1};

    self.biasLayer{1} = params{1,2};
    self.biasOutput = params{2,2};
end

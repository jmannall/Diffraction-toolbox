%% Run backwards through custome neural network to calculate derivatives and update weights and biases

function [self, averageGrad, averageSqGrad] = UpdateWeightsADAM(self, averageGrad, averageSqGrad, iteration, gradDecay, gradSqDecay, Loss)

    output = dlarray(self.outputFinal);
    if canUseGPU
        output = gpuArray(output);
    end
    [loss, delloss] = dlfeval(Loss, output, self.T);
    self.losses(end+1) = sum(loss, "all");
    
    nodeValueOutput = delloss .* self.Delactivation(self.output);
    % nodeValue is the gradient wrt the Loss at the output of the layer (before activation).
    % The nodeValue of the earlier layer is 
    nodeValueOutput = delloss;
    gradwOutput = self.hiddenOut{self.numLayers}' * nodeValueOutput;
    gradbOutput = sum(nodeValueOutput, 1);
    if self.numLayers == 1
        % nodeValueOutput * self.weightsOutput' -> the influence of the
        % output of the layer (after activation) on the next layer is
        % determined by the value of the weights (a smaller weight means a
        % smaller influence). The impact of the output of the layer (after
        % activation) on the Loss is determined by summing these weights by
        % the gradient of the node it connects to wrt to the Loss.
        nodeValueLayer{1} = nodeValueOutput * self.weightsOutput' .* self.Delactivation(self.hidden{1});
        gradwLayer{1} = self.X' * nodeValueLayer{1};
        gradbLayer{1} = sum(nodeValueLayer{1}, 1);
    else
        nodeValueLayer{self.numLayers} = nodeValueOutput * self.weightsOutput' .* self.Delactivation(self.hidden{self.numLayers});
        gradwLayer{self.numLayers} = self.hiddenOut{self.numLayers - 1}' * nodeValueLayer{self.numLayers};
        gradbLayer{self.numLayers} = sum(nodeValueLayer{self.numLayers}, 1);
        for j = 1:(self.numLayers - 2)
            i = self.numLayers - j;
            nodeValueLayer{i} = nodeValueLayer{i + 1} * self.weightsLayer{i + 1}' .* self.Delactivation(self.hidden{i});
            gradwLayer{i} = self.hiddenOut{i - 1}' * nodeValueLayer{i};
            gradbLayer{i} = sum(nodeValueLayer{i}, 1);
        end
        nodeValueLayer{1} = nodeValueLayer{2} * self.weightsLayer{2}' .* self.Delactivation(self.hidden{1});
        gradwLayer{1} = self.X' * nodeValueLayer{1};
        gradbLayer{1} = sum(nodeValueLayer{1}, 1);
    end
    params = [self.weightsLayer', self.biasLayer'; {self.weightsOutput}, {self.biasOutput}];
    grad = [gradwLayer', gradbLayer'; {gradwOutput}, {gradbOutput}];
we
    [params, averageGrad, averageSqGrad] = ADAM(params,grad,averageGrad,averageSqGrad,iteration,self.learnRate,gradDecay,gradSqDecay);
    self.weightsLayer = {params{1:self.numLayers,1}};
    self.weightsOutput = params{self.numLayers + 1,1};

    self.biasLayer = {params{1:self.numLayers,2}};
    self.biasOutput = params{self.numLayers + 1,2};
end

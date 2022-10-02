%% Run inputs through custom neural network

function self = Forward(self)
     nodeInput = self.X * self.weightsLayer{1} + self.biasLayer{1};
     self.hidden{1} = nodeInput;
     nodeOutput = self.Activation(nodeInput);
     self.hiddenOut{1} = nodeOutput;
    for i = 2:self.numLayers
         nodeInput = self.hiddenOut{i - 1} * self.weightsLayer{i} + self.biasLayer{i};
         self.hidden{i} = nodeInput;
         nodeOutput = self.Activation(nodeInput);
         self.hiddenOut{i} = nodeOutput;
    end

    self.output = nodeOutput * self.weightsOutput + self.biasOutput;
    self.outputFinal = self.Activation(self.output);
    self.outputFinal = self.output;
end
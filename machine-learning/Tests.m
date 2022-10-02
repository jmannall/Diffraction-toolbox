%% Sigmoid test
close all
x = -10:0.1:10;

y1 = Sigmoid(x);
y2 = Delsigmoid(x);

figure
plot(x, y1)
hold on
plot(x, y2)
legend('Sigmoid', 'Delsigmoid')
ylim([0 1])
grid on

%% Relu test
close all
x = -10:0.1:10;

y1 = Relu(x);
y2 = Delrelu(x);

figure
plot(x, y1)
hold on
plot(x, y2)
legend('Sigmoid', 'Delsigmoid')
ylim([-5 10])
grid on

%% Initialise
clear all

trainingData = [0, 0; 0, 1; 1, 0; 1, 1];
targetData = [0; 1; 1; 0];
learnRate = 0.2;
numEpochs = 8000;
miniBatchSize = 4;
numInput = 2;
numHidden = 2;
numOutput = 1;
range = 1e-1;
Activation = @Relu;
Delactivation = @Delrelu;

self = Initialise(trainingData, targetData, learnRate, numEpochs, miniBatchSize, numInput, numHidden, numOutput, Activation, Delactivation, range);

%% Forward
trainingData = [0, 0; 0, 1; 1, 0; 1, 1];
targetData = [0; 1; 1; 0];
learnRate = 0.2;
numEpochs = 8000;

self = struct('trainingData', trainingData, 'targetData', targetData, 'learnRate', learnRate, 'numEpochs', numEpochs, 'numLayers', 1, ...
        'weightsLayer', {{[2, -1; 2, -1]}}, 'weightsOutput', [1; 1], 'biasLayer', {{[-1, 2]}}, 'biasOutput', -1, 'losses', [], ...
        'hidden', [], 'hiddenOut', [], 'Activation', @Sigmoid, 'Delactivation', [], 'output', []);

self.T = self.targetData;
self.X = self.trainingData;

y1 = Forward(self);

%% Forward 2
trainingData = [0, 0; 0, 1; 1, 0; 1, 1];
targetData = [0; 1; 1; 0];
learnRate = 0.2;
numEpochs = 8000;
miniBatchSize = 4;
numInput = 2;
numHidden = [2, 3];
numOutput = 1;
range = 1e-1;
Activation = @Relu;
Delactivation = @Delrelu;

self = Initialise(trainingData, targetData, learnRate, numEpochs, miniBatchSize, numInput, numHidden, numOutput, Activation, Delactivation, range);

self.T = self.targetData;
self.X = self.trainingData;

y1 = Forward(self);

%% Update weights
trainingData = [0, 0; 0, 1; 1, 0; 1, 1];
targetData = [0; 1; 1; 0];
learnRate = 0.2;
numEpochs = 8000;

self = struct('trainingData', trainingData, 'targetData', targetData, 'learnRate', learnRate, 'numEpochs', numEpochs, 'numLayers', 1, ...
        'weightsLayer', {{[2, -1; 2, -1]}}, 'weightsOutput', [1; 1], 'biasLayer', {{[-1, 2]}}, 'biasOutput', -1, 'losses', [], ...
        'hidden', [], 'hiddenOut', [], 'Activation', @Sigmoid, 'Delactivation', @Delsigmoid, 'output', []);

self.T = self.targetData;
self.X = self.trainingData;

self = Forward(self);

averageGrad = [];
averageSqGrad = [];
iteration = 1;
gradDecay = 0.9;
gradSqDecay = 0.999;

[y1, averageGrad, averageSqGrad] = UpdateWeightsADAM(self, averageGrad, averageSqGrad, iteration, gradDecay, gradSqDecay, @MSELoss);

%% Update weights 2
trainingData = [0, 0; 0, 1; 1, 0; 1, 1];
targetData = [0; 1; 1; 0];
learnRate = 0.2;
numEpochs = 8000;
numInput = 2;
numHidden = [2, 3];
numOutput = 1;
range = 1e-1;
Activation = @Relu;
Delactivation = @Delrelu;

self = Initialise(trainingData, targetData, learnRate, numEpochs, numInput, numHidden, numOutput, Activation, Delactivation, range);

self.T = self.targetData;
self.X = self.trainingData;

self = Forward(self);

averageGrad = [];
averageSqGrad = [];
iteration = 1;
gradDecay = 0.9;
gradSqDecay = 0.999;

[y1, averageGrad, averageSqGrad] = UpdateWeightsADAM(self, averageGrad, averageSqGrad, iteration, gradDecay, gradSqDecay, @MSELoss);

%% Train
close all
trainingData = [0, 0; 0, 1; 1, 0; 1, 1];
targetData = [0; 1; 1; 0];
learnRate = 0.01;
numEpochs = 8000;
miniBatchSize = 4;
numInput = 2;
numHidden = 2;
numOutput = 1;
range = 1;
Activation = @Sigmoid;
Delactivation = @Delsigmoid;

self = Initialise(trainingData, targetData, learnRate, numEpochs, miniBatchSize, numInput, numHidden, numOutput, Activation, Delactivation, range);
% self.weights01 = [0.94036878, 0.10646258; 0.70646121, 0.29499702];
% self.weights12 = [0.0793157; 0.68617934];
% self.bias01 = [0.4957163, 0.22985287];
% self.bias12 = 0.796985287;
disp('Train');

%% Test
self.learnRate = 0.1;
self.outputFinal = [];
self.weights01 = self.weightsLayer{1};
self.bias01 = self.biasLayer{1};
self.weights12 = self.weightsOutput;
self.bias12 = self.biasOutput;
y1 = Train(self, @MSELoss);

x = 1:numEpochs;

figure
C = colororder;
plot(x, y1.losses, 'Color',C(2,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

%% Adam

x = {2, [2, 2]; 2, 2};
%params = {x^2, [x^3, x^2]; x^3, x};

iteration = 0;
averageGrad = [];
averageSqGrad = [];
learnRate = 0.1;
gradDecay = 0.9;
gradSqDecay = 0.99;

disp('ADAM');
for i = 1:1000
    
    grad = {2 * x{1}, [3 * x{1}, 2 * x{1}]; 3 * x{1}, 1};
    iteration = iteration + 1;
    [x, averageGrad, averageSqGrad] = ADAM(x,grad,averageGrad,averageSqGrad,iteration, learnRate, gradDecay, gradSqDecay);   
end

%% Generate Biquads
close all

numBiquads = 2;
numFilters = 256;
numFreq = 2048;

[zR, zI, pR, pI, k] = RandZPK(numBiquads, numFilters);

[tfmag, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq);

% Clip to (-128, 128) and normalise to (-1, 1)
tfmag = max(min(tfmag, 128), -128) / 128;

figure
semilogx(fvec, extractdata(tfmag))
ylim([-1 1])
xlim([20 20000])

%% Biquad Freq

n = 12;

[tfmagNBand, fvecNBand, fidx] = CreateFrequencyBands(tfmag, fvec, n);

figure
semilogx(fvec, extractdata(tfmag))
hold on
semilogx(fvecNBand, tfmagNBand, '--')

%% Train biquads
close all

trainingData = tfmag;
if canUseGPU
    trainingData = gpuArray(trainingData);
end
targetData = tfmagNBand;
learnRate = 1e-4;
numEpochs = 1000;
miniBatchSize = 1;
numInput = numFreq;
numHidden = [2048 2048];
numOutput = 4 * numBiquads + 1;

self = Initialise(trainingData, targetData, learnRate, numEpochs, miniBatchSize, numInput, numHidden, numOutput, @Relu, @Delrelu, 1e-3);
disp('Train');

y1 = Train(self, @(input, target) BiquadLoss(input, target, numBiquads, numFreq, self.miniBatchSize, fidx));

x = 1:numEpochs * self.numIterationsPerEpoch;

figure
C = colororder;
plot(x, y1.losses, 'Color',C(2,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

%% Plot

y1.X = self.trainingData;
prediction = Forward(y1);
out = prediction.outputFinal;

[zR, zI] = MinPhaseTransform(out(:,1:numBiquads)', out(:,numBiquads + 1:2 * numBiquads)');
[pR, pI] = MinPhaseTransform(out(:,2 * numBiquads + 1:3 * numBiquads)', out(:,3 * numBiquads + 1:4 * numBiquads)');
k = 100 .* Sigmoid(out(:,end))';

[tfmagPredict, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq);

figure
semilogx(fvec, 128 * extractdata(tfmag))
hold on
semilogx(fvec, extractdata(tfmagPredict), '--')
ylim([-128 128])
xlim([20 20000])

%% Create BTM responses data
close all

epsilon = 1e-3;
gParameters = struct('wedgeLength', [0.01, 50], 'wedgeIndex', [180 + epsilon, 360 - epsilon], ...
    'thetaS', epsilon, 'thetaR', epsilon, 'radiusS', [1, 1], 'radiusR', [1, 1], 'zS', 10, 'zR', 10, 'epsilon', epsilon);
numObservations = 10000;
% gParameters = struct('wedgeLength', [0.01, 20], 'wedgeIndex', [180 + epsilon, 360 - epsilon], ...
%     'thetaS', epsilon, 'thetaR', epsilon, 'radiusS', [1, 1], 'radiusR', [1, 1], 'zS', 5, 'zR', 5, 'epsilon', epsilon);
% numObservations = 50;
% fs = 48000;

[result, NNinput] = NNWedgeArray(gParameters, numObservations, fs);
tfmag = dlarray([result.tfmag])';
fvec = result(1).fvec;
tfmag = max(min(tfmag, 128), -128) / 128;

%% BTM Freq
close all

n = 6;

[tfmagNBand, fvecNBand, fidx] = CreateFrequencyBands(tfmag, fvec, n);

figure
semilogx(fvec, extractdata(tfmag))
hold on
semilogx(fvecNBand, tfmagNBand, '--')

%% Train BTM
close all

%trainingData = tfmag;
wedgeLength = max(min(NNinput.wedgeLength, 50), 0) / 50;
wedgeIndex = max(min(NNinput.wedgeIndex, 360), 180) / 180 - 1;
bendingAngle = max(min(NNinput.bendingAngle, 360), 0) / 360;
minAngle = max(min(NNinput.minAngle, 360), 0) / 360;
z1 = max(min(NNinput.z1, 50), 0) / 50;
z2 = max(min(NNinput.z2, 50), 0) / 50;
deltaZ = max(min(NNinput.deltaZ, 50), 0) / 50;
type = zeros(size(NNinput.bendingAngle));
type(NNinput.bendingAngle < 180) = 1;

trainingData = [NNinput.wedgeLength, NNinput.wedgeIndex, NNinput.bendingAngle, NNinput.minAngle, NNinput.z1, NNinput.z2, NNinput.deltaZ, type];
if canUseGPU
    trainingData = gpuArray;
end
targetData = tfmagNBand;
learnRate = 1e-2;
numEpochs = 500; 
miniBatchSize = 128;
numFreq = size(fidx, 2);
numInput = size(trainingData, 2);
numHidden = [8 8];
numBiquads = 2;
numOutput = 4 * numBiquads + 1;

self = Initialise(trainingData, targetData, learnRate, numEpochs, miniBatchSize, numInput, numHidden, numOutput, @Relu, @Delrelu, 1e-1);
disp('Train');

y1 = Train(self, @(input, target) BiquadLoss(input, target, numBiquads, numFreq, self.miniBatchSize, fidx));

x = 1:length(y1.losses);

figure
C = colororder;
plot(x, y1.losses, 'Color',C(2,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

%% Plot

y1.X = self.trainingData;
y1.T = self.targetData;
prediction = Forward(y1);
out = prediction.outputFinal;

[zR, zI] = MinPhaseTransform(out(:,1:numBiquads)', out(:,numBiquads + 1:2 * numBiquads)');
[pR, pI] = MinPhaseTransform(out(:,2 * numBiquads + 1:3 * numBiquads)', out(:,3 * numBiquads + 1:4 * numBiquads)');
k = 100 * Sigmoid(out(:,end))';

[tfmagPredict, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq);
test = tfmag(1:100,:);
figure
semilogx(fvec, 128 * extractdata(tfmag(150:200,:)))
hold on
semilogx(fvec, extractdata(tfmagPredict(150:200,:)), '--')
ylim([-128 0])
xlim([20 20000])

%% Create BTM responses data (testing)

epsilon = 1e-3;
gParameters = struct('wedgeLength', [0.01, 20], 'wedgeIndex', [180 + epsilon, 360 - epsilon], ...
    'thetaS', [epsilon, 999], 'thetaR', [epsilon, 999], 'radiusS', [1, 1], 'radiusR', [1, 1], 'zS', [-5, 999], 'zR', [-5, 999], 'epsilon', epsilon);
numObservations = 40;
fs = 48000;

[geometry, NNinput] = NNGeometry(gParameters, numObservations);
[result, NNinput] = NNWedgeArray(geometry, NNinput, gParameters, numObservations, fs);
tfmag = dlarray([result.tfmag])';
fvec = result(1).fvec;
tfmag = max(min(tfmag, 128), -128) / 128;

%% BTM Freq
close all

n = 12;

[tfmagNBand, fvecNBand, fidx] = CreateFrequencyBands(tfmag, fvec, n);

for i = 1:numObservations
    figure(1)
    semilogx(fvec, extractdata(tfmag(i,:)))
    hold on
    semilogx(fvecNBand, tfmagNBand(i,:), '--')
    hold off
    x = [i, NNinput.z1(i), NNinput.z2(i), NNinput.wedgeLength(i)]
end

%% Test BTM Network plot

trainingData = tfmag;
if canUseGPU
    trainingData = gpuArray;
end
targetData = tfmagNBand;

y1.X = trainingData;
y1.T = targetData;
tests = Forward(y1);
out = tests.outputFinal;

[zR, zI] = MinPhaseTransform(out(:,1:numBiquads)', out(:,numBiquads + 1:2 * numBiquads)');
[pR, pI] = MinPhaseTransform(out(:,2 * numBiquads + 1:3 * numBiquads)', out(:,3 * numBiquads + 1:4 * numBiquads)');
k = 100 * Sigmoid(out(:,end))';

[tfmagPredict, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq);

figure
semilogx(fvec, 128 * extractdata(tfmag))
hold on
semilogx(fvec, extractdata(tfmagPredict), '--')
ylim([-128 0])
xlim([20 20000])
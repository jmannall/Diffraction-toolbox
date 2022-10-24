clc;
clear;
close all;
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultLineMarkerSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NN

% Hyperparamters
% learnRate = 1e-2;
% gradDecay = 0.9;
% sqGradDecay = 0.999;
% maxGradient = 5;

% x.lR = learnRate;
% x.gD = gradDecay;
% x.sGD = sqGradDecay;
% x.mG = maxGradient;
% loss = XORNeuralNetwork(x);

learnRate = optimizableVariable('lR',[1e-5 1e-1],'Type','real','Transform','log');
gradDecay = optimizableVariable('gD',[0.5 1],'Type','real');
sqGradDecay = optimizableVariable('sGD',[0.8 1],'Type','real');
maxGradient = optimizableVariable('mG',[0.5 10],'Type','real','Transform','log');

func = @(x)XORNeuralNetwork(x);
result = bayesopt(func, [learnRate, gradDecay, sqGradDecay, maxGradient], 'UseParallel',true);

xObs = result.XAtMinObjective;
xEst = result.XAtMinEstimatedObjective;

[lossObs, netObs] = XORNeuralNetwork(xObs);
[lossEst, netEst] = XORNeuralNetwork(xEst);

%% XOR Network
    trainingData = ([0 0
    1 0
    0 1
    1 1])';

targetData = [1 0 0 1];

[numInputs, dataSize] = size(trainingData);
numLayers = 1;
hiddenLayerSize = 2;
numOutputs = size(targetData, 1);
alpha = 0.2;

net = CreateNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha);
learnables = net.Learnables.Value;

miniBatchSize = 4;
numEpochs = 500;
numIterationsPerEpoch = floor(dataSize./miniBatchSize);

[lineIterationLoss, lineEpochLoss] = CreateAnimatedLinePlot();

averageGrad = [];
averageSqGrad = [];

iteration = 0;
start = tic;

lossFunc = @(net, trainingData, targetData) MeanSquaredError(net, trainingData, targetData);
losses = zeros(1, numEpochs);

tic
for epoch = 1:numEpochs
    % Shuffle data.
    idx = randperm(size(trainingData, 2));
    targetData = targetData(:,idx);
    trainingData = trainingData(:,idx);
    for i = 1:numIterationsPerEpoch
        iteration = iteration + 1;

        [X, T] = ProcessNNBatchInputData(trainingData, targetData, miniBatchSize, i);

        % Evaluate the model loss and gradients using dlfeval and the
        % modelLoss function.input
        [loss, state, gradients] = dlfeval(lossFunc,net,X,T);
        % Update normalisation paramters
        net.State = state;

        % Update the network parameters using the Adam optimizer.
        % Clip gradients
        gradients = ClipGradients(gradients, maxGradient);
        [net,averageGrad,averageSqGrad] = adamupdate(net,gradients,averageGrad,averageSqGrad,iteration,learnRate,gradDecay,sqGradDecay);
        
        idx = numIterationsPerEpoch * (epoch - 1) + i;
        loss = double(loss);
        losses(idx) = loss;
    end

    UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)
end
toc



%% Generate BTM Data
close all
disp('Create BTM data');

% fs = 48000;
% 
% epsilon = 1e-3;
% gParameters = struct('wedgeLength', 20, 'wedgeIndex', [180 + epsilon, 360 - epsilon], ...
%     'thetaS', [epsilon, 999], 'thetaR', [epsilon, 999], 'radiusS', 1, 'radiusR', 1, 'zS', 10, 'zR', 10, 'epsilon', epsilon);
% numObservations = 20000;
% 
% [geometry, NNinput] = NNGeometry(gParameters, numObservations);
% [result, NNinput] = NNWedgeArray(geometry, NNinput, gParameters, numObservations, fs);

numObservations = 274;
wedgeIndex = 359;
bendingAngle = 5.5:355.5;
minAngle = 40.01;

geometry = Geometry(wedgeIndex, bendingAngle, minAngle);

disp('Run')
[result, geometry] = SingleWedgeArray(geometry, 20,1,1,10,10,48000);
bendingAngle = geometry.bendingAngle;

tfmag = dlarray([result.tfmag])';
tfmag = max(min(tfmag, 128), -128) / 128;
fvec = result(1).fvec;
tfmagBTM = extractdata(tfmag);
tfcomplex = [result.tfcomplex]';

numBiquads = 2;
numFreq = size(tfmag, 2);

%% Biquad Freq

n = 12;

[tfmagNBand, fvecNBand, fidx] = CreateFrequencyBands(tfmag, fvec, n);

figure
semilogx(fvec, extractdata(tfmag(1:30:270,:)))
hold on
semilogx(fvecNBand, tfmagNBand(1:30:270,:), '--')

%% Initialise network
disp('Intialise network');
close all

% percent = 0.9;
% split = numObservations * percent;
% idx = 1:split;
% trainingData = [NNinput.wedgeIndex(idx,:), NNinput.bendingAngle(idx,:), NNinput.minAngle(idx,:)]';
% targetData = struct('mag', tfmagNBand(idx,:)', 'phase', angle(tfcomplex(idx,:)'));
% targetBTM = tfmagBTM(idx,:)';

const = ones(numObservations, 1);
wedgeI = const * wedgeIndex;
minA = const * minAngle;
bendingA = bendingAngle;
trainingData = [wedgeI, bendingA, minA]';
targetData = tfmagNBand';
phaseData = dlarray(angle(tfcomplex)');
targetBTM = tfmagBTM;

% idx = split + 1;
% inputData = [NNinput.wedgeIndex(idx:end), NNinput.bendingAngle(idx:end), NNinput.minAngle(idx:end)]';
% testData = tfmagNBand(idx:end,:)';
% testBTM = tfmagBTM(idx:end,:)';

inputData = [wedgeI, bendingA, minA]';
testData = tfmagNBand';
testBTM = tfmagBTM;
testDataPhase = phaseData;

numInputs = size(trainingData, 1);
numOutputs = 4 * numBiquads + 1 + 1;
hiddenLayerSize = 32;
alpha = 0.2;

lgraph = layerGraph();

tempLayers = [
    featureInputLayer(numInputs,"Name","InputLayer")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerOne")
    batchNormalizationLayer("Name", "NormalisationLayerOne")
    leakyReluLayer(alpha, "Name", "ActivationLayerOne")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerTwo")
    batchNormalizationLayer("Name", "NormalisationLayerTwo")
    leakyReluLayer(alpha, "Name", "ActivationLayerTwo")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerThree")
    batchNormalizationLayer("Name", "NormalisationLayerThree")
    leakyReluLayer(alpha, "Name", "ActivationLayerThree")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerFour")
    batchNormalizationLayer("Name", "NormalisationLayerFour")
    leakyReluLayer(alpha, "Name", "ActivationLayerFour")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerFive")
    batchNormalizationLayer("Name", "NormalisationLayerFive")
    leakyReluLayer(alpha, "Name", "ActivationLayerFive")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerSix")
    batchNormalizationLayer("Name", "NormalisationLayerSix")
    leakyReluLayer(alpha, "Name", "ActivationLayerSix")
    fullyConnectedLayer(numOutputs,"Name","OutputLayer")
];

lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

plot(lgraph);

net = dlnetwork(lgraph);

learnables = net.Learnables.Value;

percent = 1;
miniBatchSize = 128;
numEpochs = 100;
numIterationsPerEpoch = floor(percent * numObservations./miniBatchSize);
figure
C = colororder;
lineLossTrain = animatedline('Color',C(2,:));
lineLossEpoch = animatedline('Color',C(1,:));
lineLossTest = animatedline('Color',C(3,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

averageGrad = [];
averageSqGrad = [];

iteration = 0;
start = tic;

fs = 48000;

lossFunc = @(net, inputData, targetData, phaseData) NNBiquadLossPhase(net, inputData, targetData, phaseData, numBiquads, numFreq, miniBatchSize, fidx);
lossFunc_Test = @(net, inputData, targetData) NNBiquadLoss_Test(net, inputData, targetData, numBiquads, numFreq, miniBatchSize, fidx);
losses = zeros(1,numIterationsPerEpoch * numEpochs);
[epochLosses, testLosses] = deal(zeros(1,numEpochs));

tic
for epoch = 1:numEpochs
    % Shuffle data.
    idx = randperm(size(trainingData, 2));
    targetData = targetData(:,idx);
    phaseData = phaseData(:,idx);
    trainingData = trainingData(:,idx);
    targetBTM = targetBTM(:,idx);
    for i = 1:numIterationsPerEpoch
        iteration = iteration + 1;

        % Read mini-batch of data and convert the labels to dummy
        % variables.
        idx = (i-1)*miniBatchSize+1:i*miniBatchSize;
        X = trainingData(:,idx);
        T = targetData(:,idx);
        P = phaseData(:,idx);

%         T = zeros(numOutputs, miniBatchSize,"single");
%         for c = 1:numClasses
%             T(c,inputTraingData(idx)==classes(c)) = 1;
%         end

        % Convert mini-batch of data to a dlarray.
        X = dlarray(single(X), "CB");

        % If training on a GPU, then convert data to a gpuArray.
%         if canUseGPU
%             executionEnvironment = "gpu";
%             numberOfGPUs = gpuDeviceCount("available");
%             pool = parpool(numberOfGPUs);
%             X = gpuArray(X);
%         end

        % Evaluate the model loss and gradients using dlfeval and the
        % modelLoss function.input
        [loss, state, gradients] = dlfeval(lossFunc,net,X,T,P);
        net.State = state;
        %loss = lossFunc_Test(net,X,T);

        t = 80;
        alpha = min(iteration / t, 1);
        learnRateStart = 1e-2;
        learnRateEnd = 1e-4;
        learnRate = (1 - alpha) * learnRateStart + alpha * learnRateEnd;
        %learnRate = 1e-3;
        gradDecay = 0.9;
        sqGradDecay = 0.999;
        % Update the network parameters using the Adam optimizer.
        gradValue = gradients.Value;
        maxGradient = 0.9;
        for j = 1:length(gradValue)
            clip = false;
            gradient = gradValue{j};
            for k = 1:size(gradient, 1)
                g = extractdata(gradient(k,:));
                gNorm = norm(g);
                if gNorm > maxGradient
                    clip = true;
                    gradient(k,:) = g * maxGradient ./ gNorm;
                end
            end
            if clip
                gradValue{j} = gradient;
            end
        end
        gradients.Value = gradValue;
        [net,averageGrad,averageSqGrad] = adamupdate(net,gradients,averageGrad,averageSqGrad,iteration,learnRate,gradDecay,sqGradDecay);
        
        idx = numIterationsPerEpoch * (epoch - 1) + i;
        loss = double(loss);
        losses(idx) = loss;
        if mod(iteration, 100) == 0
            D = duration(0,0,toc(start),'Format',"hh:mm:ss");
            addpoints(lineLossTrain,idx,loss)
            title("Epoch: " + epoch + ", Elapsed: " + string(D))
            drawnow
            disp(['Loss: ', num2str(loss)]);
        end
    end
    [testLoss, testLossAll] = NNBiquadLossPhase_Test(net, dlarray(single(inputData), "CB"), testData, testDataPhase, numBiquads, numFreq, numObservations, fidx);
    testLosses(epoch) = testLoss;
    idx = randperm(size(trainingData, 2));
    idx = idx(1:round(numObservations));
    [epochLoss, epochLossAll] = NNBiquadLossPhase_Test(net, dlarray(single(trainingData(:,idx)), "CB"), targetData(:,idx), phaseData(:,idx), numBiquads, numFreq, numObservations, fidx);
    epochLosses(epoch) = epochLoss;
    idx = numIterationsPerEpoch * epoch;
    D = duration(0,0,toc(start),'Format',"hh:mm:ss");
    addpoints(lineLossEpoch,idx,epochLoss)
    addpoints(lineLossTest,idx,testLoss)
    title("Epoch: " + epoch + ", Elapsed: " + string(D))
    drawnow
    disp([num2str(epoch), ' Epoch loss: ', num2str(epochLoss)]);
    disp(['Test loss: ', num2str(testLoss)]);
    avEpochLoss = sum(epochLosses((max(1,epoch - 9)):epoch)) / 10;
end
toc

epochLossAll = extractdata(epochLossAll);
testLossAll = extractdata(testLossAll);

%% Save
disp('Save network')

save('NeuralNetwork_TestingP003.mat', 'net', 'losses', 'epochLosses', 'testLosses', 'epochLossAll', 'testLossAll', '-v7.3')

%% Plot

disp('Plot training losses')

xObs = 1:numEpochs * numIterationsPerEpoch;
y = 1:numIterationsPerEpoch:numEpochs * numIterationsPerEpoch;

figure
C = colororder;
plot(xObs, losses, 'Color',C(2,:));
hold on
plot(y, epochLosses, 'Color',C(1,:));
plot(y, testLosses, 'Color',C(3,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

disp('Training Complete');

%% Analyse losses
close all
disp('Analyse Losses')

load('NeuralNetwork_014.mat', 'net', 'epochLosses', 'testLosses', 'epochLossAll', 'testLossAll');

figure
histogram(epochLossAll)
xlim([0 30])
ylim([0 1000])
title('Training data')
xlabel('Frequency')
ylabel('Root mean squared error')

meanEpochLoss = mean(epochLossAll);
medianEpochLoss = median(epochLossAll);
maxEpochLoss = max(epochLossAll);
minEpochLoss = min(epochLossAll);

figure
histogram(testLossAll)
xlim([0 30])
ylim([0 1000])
title('Testing data')
xlabel('Frequency')
ylabel('Root mean squared error')

meanTestLoss = mean(testLossAll);
medianTestLoss = median(testLossAll);
maxTestLoss = max(testLossAll);
minTestLoss = min(testLossAll);

%% Test network
disp('Test Network');

output = predict(net, dlarray(inputData, "CB"));

[zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);
[tfmag, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq);

tfmag = extractdata(tfmag)';

%% Plot
for i = 1:numObservations * 0.1
    figure(4)
    semilogx(fvec, tfmag(:,i))
    hold on
    semilogx(fvec, 128 * testBTM(:,i))
    semilogx(fvecNBand, 128 * testData(:,i),'--')
    hold off
    xlim([20, 20000])
    ylim([-60 0])
    legend('Prediction', 'Target')
end
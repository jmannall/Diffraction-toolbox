%% Train a neural network to match IIR filters to input frequency responses.
%% Based on paper " Direct Design of Biquad Filter Cascades with Deep
%% Learning by Sampling Random Polynomials" by Joseph T Colonel.

clc;
clear;
close all;
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultLineMarkerSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Initialise network
disp('Intialise network');
close all

numInputs = numFreq;
numOutputs = 4 * numBiquads + 1;
hiddenLayerSize = 2048;
alpha = 0.2;

lgraph = layerGraph();

tempLayers = [
    featureInputLayer(numInputs,"Name","InputLayer")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerOne")
    leakyReluLayer(alpha, "Name", "ActivationLayerOne")
    fullyConnectedLayer(hiddenLayerSize,"Name","HiddenLayerTwo")
    leakyReluLayer(alpha, "Name", "ActivationLayerTwo")
    fullyConnectedLayer(numOutputs,"Name","OutputLayer")
];

lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

plot(lgraph);

net = dlnetwork(lgraph);

learnables = net.Learnables.Value;

trainingData = tfmag';
targetData = tfmagNBand';
miniBatchSize = 128;
numEpochs = 1000;
numObservations = numFilters;
numIterationsPerEpoch = floor(numObservations./miniBatchSize);

figure
C = colororder;
lineLossTrain = animatedline('Color',C(2,:));
lineLossEpoch = animatedline('Color',C(1,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

averageGrad = [];
averageSqGrad = [];

iteration = 0;
start = tic;

fs = 48000;

lossFunc = @(net, inputData, targetData) NNBiquadLoss(net, inputData, targetData, numBiquads, numFreq, miniBatchSize, fidx);
losses = zeros(1,numIterationsPerEpoch * numEpochs);
epochLosses = zeros(1,numEpochs);

for epoch = 1:numEpochs
    % Shuffle data.
    idx = randperm(size(trainingData, 2));
    targetData = targetData(:,idx);
    trainingData = trainingData(:,idx);

    for i = 1:numIterationsPerEpoch
        iteration = iteration + 1;

        % Read mini-batch of data and convert the labels to dummy
        % variables.
        idx = (i-1)*miniBatchSize+1:i*miniBatchSize;
        X = trainingData(:,idx);
        T = targetData(:,idx);

        % Convert mini-batch of data to a dlarray.
        X = dlarray(single(X), "CB");

        % If training on a GPU, then convert data to a gpuArray.
        if canUseGPU
            X = gpuArray(X);
        end

        % Evaluate the model loss and gradients using dlfeval and the
        % modelLoss function.input
        [loss,gradients] = dlfeval(lossFunc,net,X,T);

        learnRate = 1e-4;
        gradDecay = 0.9;
        sqGradDecay = 0.999;
        % Update the network parameters using the Adam optimizer.
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
    epochLoss = sum(losses(((epoch - 1) * numIterationsPerEpoch + 1):epoch * numIterationsPerEpoch)) / numIterationsPerEpoch;
    idx = numIterationsPerEpoch * epoch;
    epochLosses(epoch) = epochLoss;
    D = duration(0,0,toc(start),'Format',"hh:mm:ss");
    addpoints(lineLossEpoch,idx,epochLoss)
    title("Epoch: " + epoch + ", Elapsed: " + string(D))
    drawnow
    disp(['Epoch loss: ', num2str(epochLoss)]);
    avEpochLoss = sum(epochLosses((max(1,epoch - 9)):epoch)) / 10;
end


%% Plot

x = 1:numEpochs * numIterationsPerEpoch;
y = 1:numIterationsPerEpoch:numEpochs * numIterationsPerEpoch;

figure
C = colororder;
plot(x, losses, 'Color',C(2,:));
hold on
plot(y, epochLosses, 'Color',C(1,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

disp('Training Complete');

%% Test network
disp('Test Network');

output = predict(net, dlarray(trainingData, "CB"));

[zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);
[tfmag, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq);

tfmag = extractdata(tfmag)';

close all
for i = 1:size(tfmag,2)
    
    figure(3)
    semilogx(fvec, tfmag(:,i))
    hold on
    semilogx(fvec, 128 * extractdata(trainingData(:,i)))
    semilogx(fvecNBand, 128 * targetData(:,i),'--')
    hold off
    xlim([20, 20000])
    ylim([-128 128])
    legend('Prediction', 'Target')
end
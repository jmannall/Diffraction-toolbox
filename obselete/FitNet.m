numData = 10000;

data1 = randi([0 1],[1 numData]);
data2 = randi([0 1],[1 numData]);

XORresult = xor(data1, data2);

inputData = [data1; data2];
targetData = double(XORresult);

inputData = dlarray(inputData,"CB");

numInputs = size(inputData, 1);
numOutputs = 1;

hiddenLayerSize = 2;

net = fitnet(hiddenLayerSize);

net.performFcn = 'mymse';
view(net)
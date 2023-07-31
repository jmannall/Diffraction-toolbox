function [z, p, k] = mySmallNN(in)

    persistent mynet;
    if isempty(mynet)
        mynet = coder.loadDeepLearningNetwork('NNSaves_FinalRun/Run6/IIR-4_20_0001.mat');
    end

    in = dlarray(single(in'), "CB");
    out = extractdata(predict(mynet,in));
    [z, p, k] = TransformNNOutput(out);
end
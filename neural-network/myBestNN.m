function [z, p, k] = myBestNN(in)

    persistent mynet;
    if isempty(mynet)
        mynet = coder.loadDeepLearningNetwork('NNSaves_FinalRun/Run6/IIR-7_36_0001.mat');
    end
    
    in = dlarray(single(in'), "CB");
    out = extractdata(predict(mynet,in));
    [z, p, k] = TransformNNOutput(out);
end
function [loss, gradients] = TransformationLoss(net, psi, psiInv, x, y, NNlossFunc, numIIRFilters, gamma)

    NNLoss = NNlossFunc(net, x, y);
    w = psi * stripdims(x);
    xHat = psiInv * w;
    r = x - xHat;
    input = xHat + extractdata(r);

    input = dlarray(input, "CB");
    output = predict(net, input);
    [z, p, k] = CreateIIRFromNNOutput(output, numIIRFilters);

    z = mean(z, 2);
    p = mean(p, 2);
    k = mean(k);
   
%     TRIM = dlarray (zeros(1));
%     input = num2cell([z; p; k]);
%     for i = 1:size(input, 1)
%         attrib = dlgradient(input{i}, w);
%         TRIM = TRIM + mean(abs(attrib), 'all');
%     end
%     TRIM = TRIM / size(input, 1);
    TRIM = dlarray (zeros(1));
    input = num2cell([z; p; k]);
    for i = 1:length(input)
        attrib = dlgradient(input{i}, w, 'EnableHigherDerivatives',true);
        TRIM = TRIM + sum(abs(attrib), 'all');
    end
    loss = NNLoss + gamma * TRIM;
    gradients = dlgradient(loss, psi);
end
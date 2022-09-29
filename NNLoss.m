function [loss, gradients] = NNLoss(net, input, tfvalue)

    numFreq = size(tfvalue,1);
    numObservations = size(tfvalue,2);

    % Get outputs of the NN
    [zR, zI, pR, pI, k] = forward(net,input);
    
    k = k * 0.1;
    [tfmag, fvec] = CreateBiquad(zR, zI, pR, pI, k);

    loss = sum((tfmag - tfvalue) .^ 2,"all") / (numFreq * numObservations);

%     f4 = figure(4);
%     clf(f4)
%     semilogx(fvec, tfvalue)
%     hold on
%     semilogx(fvec, extractdata(tfmag), '--')
%     hold off

    disp(['Loss: ', num2str(loss)]);
    gradients = dlgradient(loss,net.Learnables);

end
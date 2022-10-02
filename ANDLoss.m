function [loss, gradients, accuracy] = ANDLoss(net, input, target)
       
    % Get outputs of the NN
    output = forward(net,input);
    loss = mse(output, target);
    gradients = dlgradient(loss,net.Learnables);
    
    accuracy = mean(round(output)==target);

    disp(['Loss: ', num2str(loss)]);

end
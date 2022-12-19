function cost = CalculateNNIIRCost(nL, nN, nI, nO, gx)

    % gx = 5; % Activation layer and normalisation cost per node
    % nN = 42; % Number of nodes per layer
    % nL = 5; % Number of layers
    % nO = 5; % Number of outputs
    % nI = 8; % Number of inputs
    CcPZTransform = 5;
    numPZ = nO - 1;
    CcKTransform = 1;
    CcOutputTransform = numPZ * CcPZTransform + CcKTransform;
    
    nNTotal = nN * nL + nO; % Total number of nodes in DNN (i.e the number of bias additions)
    nNActivations = nN * nL; % Total number of nodes in DNN that require activation/normalisation
    nCTotal = nI * nN + nN ^ 2 * (nL - 1) + nN * nO; % Total number of connections in DNN (i.e the number of weights)
    
    CcNN = nNTotal + nNActivations * gx + nCTotal + CcOutputTransform;
    CcIIR = 6;
    
    cost = CcNN + CcIIR;
end
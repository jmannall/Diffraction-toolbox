function H = SingleHpm(z, f, gParameters, controlparameters) % 161

    thetaPlus = gParameters.thetaS + gParameters.thetaR; % 1
    thetaMinus = gParameters.thetaR - gParameters.thetaS; % 1
    H = SingleH(z, thetaPlus, f, gParameters, controlparameters) + SingleH(z, thetaMinus, f, gParameters, controlparameters); %2n + 1 -> 159
end
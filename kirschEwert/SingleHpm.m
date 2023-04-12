function H = SingleHpm(z, f, gParameters, controlparameters)

    thetaPlus = gParameters.thetaS + gParameters.thetaR;
    thetaMinus = gParameters.thetaR - gParameters.thetaS;
    H = SingleH(z, thetaPlus, f, gParameters, controlparameters) + SingleH(z, thetaMinus, f, gParameters, controlparameters);
end
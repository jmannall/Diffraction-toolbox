function output = CreateNNOutput(net, wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR)

    X = NNInputFromGeometry(wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR);
    
    output = predict(net, X);
end
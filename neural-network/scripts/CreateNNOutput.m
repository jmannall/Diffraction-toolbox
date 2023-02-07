function output = CreateNNOutput(net, wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR, doReflection)

    X = NNInputFromGeometry(wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR, doReflection);
    
    output = predict(net, X);
end
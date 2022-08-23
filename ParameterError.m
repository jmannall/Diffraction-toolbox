function z = ParameterError(x, tfvalue, geometry, fs)

    numResults = length(BTMresult);
    w = geometry.wedge;
    mA = geometry.minAngle;
    bA = geometry.bendingAngle;
    
    z = 0;
    
    for i = 1:numResults
        pole1 = f(w, mA, mB);
        pole2 = g(w, mA, mB);
        zero1 = h(w, mA, mB);
        zero2 = x(w, mA, mB);
        gain = y(w, mA, mB);
    
        p = [pole1; pole2];
        z = [zero1; zero2];
        k = gain;
    
        [tf, fvec] = IIRFilter(p, z, k, fs);
        z = z + Error(tf, tfvalue, fvec);
    end
    
    z = z / numResults;
end
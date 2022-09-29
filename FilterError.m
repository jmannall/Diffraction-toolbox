function z = FilterError(x, tfvalue, fs)

    % Create zeros, poles and gain from x   
    z = [x(1); x(3)];
    p = [x(2); x(4)];
    k = x(5);

    % Create tf and fvec of the IIR filter
    [tf, fvec] = IIRFilter(z, p, k, fs);

    % Find the MSE compared to BTM response tfvalue
    z = Error(tf, tfvalue, fvec);
end

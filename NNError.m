function MSE = NNError(z, p, k, tfvalue, fs)

    numGeom = length(z);

    % Create tf and fvec of the IIR filters
    MSE = 0;
    for i = 1:numGeom
        [tfmag, fvec] = NNIIRFilter(z(:,i), p(:,i), k(i), fs);

        % Find the MSE compared to BTM response tfvalue
        MSE = MSE + Error(tfmag, tfvalue(:,i), fvec);
    end
    MSE = MSE / numGeom;
end
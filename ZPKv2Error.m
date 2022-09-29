function MSE = ZPKv2Error(x, input, tfvalue, fs)
   
    numGeom = length(input);

    prediction = x * input;

    z = [prediction(1,:); prediction(3,:)];
    p = [prediction(2,:); prediction(4,:)];
    k = prediction(5,:);

    % Create tf and fvec of the IIR filters
    MSE = 0;
    for i = 1:numGeom

        [tf, fvec] = IIRFilter(z(:,i), p(:,i), k(i), fs);

        % Find the MSE compared to BTM response tfvalue
        MSE = MSE + Error(tf, tfvalue(:,i), fvec);
    end

    MSE = MSE / numGeom;
end
function x = ZPKConstraints(x, input)
    
    prediction = (input*x)';

    indexMin = prediction <= -1;

    filterParamterMin = max(min(sum(indexMin, 2), 1), 0);

    index = find(filterParamterMin == 1);
    for i = 1:sum(filterParamterMin)
        [~, location] = max(prediction(index(i),:));
        w = input(location, 2);
        bA = input(location, 3);
        mA = input(location, 4);
        x(1,index(i)) = max(x(1,index(i)), -1);
        scale = (-0.99 - x(1,index(i))) / (w * x(2,index(i)) + bA * x(3,index(i)) + mA * x(4,index(i)));
        x(2:4,index(i)) = x(2:4,index(i)) * scale;
    end

    prediction = (input*x)';

    indexMax = prediction >= 1;
    filterParamterMax = max(min(sum(indexMax, 2), 1), 0);
    index = find(filterParamterMax == 1);
    for i = 1:sum(filterParamterMax)
        [~, location] = max(prediction(index(i),:));
        w = input(location, 2);
        bA = input(location, 3);
        mA = input(location, 4);
        x(1,index(i)) = max(x(1,index(i)), -1);
        scale = (0.99 - x(1,index(i))) / (w * x(2,index(i)) + bA * x(3,index(i)) + mA * x(4,index(i)));
        x(2:4,index(i)) = scale * x(2:4,index(i));
    end

    prediction = (input*x)';
end
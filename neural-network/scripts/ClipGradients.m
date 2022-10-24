function gradients = ClipGradients(gradients, v)

    gradValue = gradients.Value;

    for j = 1:length(gradValue)
        clip = false;
        gradient = gradValue{j};
        for k = 1:size(gradient, 1)
            g = extractdata(gradient(k,:));
            gNorm = norm(g);
            if gNorm > v
                clip = true;
                gradient(k,:) = g * v ./ gNorm;
            end
        end
        if clip
            gradValue{j} = gradient;
        end
    end
    gradients.Value = gradValue;
end
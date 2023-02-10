function vector = NormaliseVector(vector)
    vector = vector ./ vecnorm(vector, 2, 2);
end
%% Find the angle in degrees between two vectors

function theta = AngleBetweenVectors(a, b)
    theta = acosd(dot(a, b) ./ (vecnorm(a) .* vecnorm(b)));
end
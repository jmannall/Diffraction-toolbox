function theta = AngleBetweenVectors(a, b)
    theta = acosd(dot(a, b) ./ (vecnorm(a) .* vecnorm(b)));
end
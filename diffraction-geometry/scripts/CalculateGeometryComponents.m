function [radius, theta, z] = CalculateGeometryComponents(position, wedgeIndex)

    radius = vecnorm(position(:,1:2), 2, 2);
    theta = wedgeIndex + sign(position(:,2)) .* acosd(position(:,1) ./ radius) - 180;
    z = position(:,3);
end
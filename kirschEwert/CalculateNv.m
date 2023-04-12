function Nv = CalculateNv(v, theta)
    Nv = (v * sqrt(1 - cos(v * pi) * cos(v * theta))) / (cos(v * pi) - cos(v * theta));
end
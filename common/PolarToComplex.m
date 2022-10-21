%% Convert polar coordinates to a complex number

function z = PolarToComplex(r, arg)
    a = r.*cos(arg);
    b = r.*sin(arg);

    z = a+b*1i;
end
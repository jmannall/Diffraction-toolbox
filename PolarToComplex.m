function output = PolarToComplex(r, arg)
    a = r*cos(arg);
    b = r*sin(arg);

    output = a+b*1i;
end
function [r, arg] = ComplexToPolar(z)
    r = abs(z);
    arg = angle(z);
end
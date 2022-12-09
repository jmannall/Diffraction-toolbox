function [b, a] = IIRFilterParameterCoefficients(lpFc, hsFc, G, k, fs)

    [b, a] = deal(zeros(2, 2, length(lpFc)));
    [b(:,1,:), a(:,1,:)] = LowPassCoefficients(lpFc, fs, k);
    [b(:,2,:), a(:,2,:)] = HighShelfCoefficients(hsFc, G, fs);

    b = real(b);
    a = real(a);
    if isDlarray
        b = stripdims(b);
        a = stripdims(a);
    end
end
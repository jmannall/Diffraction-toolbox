function [b, a] = IIRFilterParameterCoefficients(lpFc, hsFc, G, k, fs)

    isDlarray = isdlarray(lpFc);

    if isDlarray
        [b, a] = deal(dlarray(zeros(2, 2, length(lpFc))));
    else
        [b, a] = deal(zeros(2, 2, length(lpFc)));
    end
    [b(:,1,:), a(:,1,:)] = LowPassCoefficients(lpFc, fs, k);
    [b(:,2,:), a(:,2,:)] = HighShelfCoefficients(hsFc, G, fs);
    b = real(b);
    a = real(a);
    if isDlarray
        b = stripdims(b);
        a = stripdims(a);
    end
end
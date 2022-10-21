function output = EquationQuarter(beta, n, k, L, plus, B) % 35

    cotArg = (pi + PlusMinus(beta, plus)) / (2 * n); % 4
    if (abs(cotArg) < 0.01 ) % Going to be singular
        disp('Almost singular mode');
        [~, betaarg] = Apm(n, beta, plus);
        eps = pi +  betaarg;
        if eps == 0
            eps = 0.001;
        end
        if (abs(eps) >= 1)
            msg = 'eps more than 1';
            error(msg)
        end
        output = n * exp(1i * pi / 4) * ((sqrt(2 * pi * k * L) * sign(eps) - 2 * k * L * eps * exp(1i * pi / 4)));
    else
        output = cot(cotArg) * FuncF(B * k * L * Apm(n, beta, plus)); % 2 + FuncF + Apm -> 31
    end
end

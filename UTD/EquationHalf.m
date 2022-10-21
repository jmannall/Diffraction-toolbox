function output = EquationHalf(beta, n, k, L)

    output = EquationQuarter(beta, n, k, L, true) + EquationQuarter(beta, n, k, L, false);
%     output = cot((pi + beta) / (2 * n)) .* FuncF(k * L * Apm(n, beta, true)) + ...
%     cot((pi - beta) / (2 * n)) .* FuncF(k * L * Apm(n, beta, false));
end
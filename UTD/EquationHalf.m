function output = EquationHalf(beta, n, k, L, B) % 71

    output = EquationQuarter(beta, n, k, L, true, B) + EquationQuarter(beta, n, k, L, false, B); % 1 + 2 * eqQuarter -> 71
end
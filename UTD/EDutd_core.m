%% UTD model

function D = EDutd_core(thetai, li, thetad, ld, v, phii, k, softboundary)
    A = EDkouyoumjianlenatten(li, ld, k);
    frontfactor = -A * v * exp(-1i * pi / 4) / (2 * sqrt(2 * pi * k) * sin(phii));
    Lambda = li * ld * (sin(phii) ^ 2) / (li + ld);
    n = 1 / v;  
    function ret = equationhalf(beta)
        function ret = equationquarter(beta, plus)
            function ret = Npm(plus, beta)
                if (plus)
                    if (beta <= pi * (n - 1))
                        ret = 0;
                    else
                        ret = 1;
                    end
                else
                    if (beta < pi * (1 - n))
                        ret = -1;
                    elseif ( beta > pi * (1 + n))
                        ret = 1;
                    else
                        ret = 0;
                    end
                end
            end
            betaarg = beta - (2 * pi * n * Npm(plus, beta));
            cotarg = 0.5 * v * (pi + EDpm(plus, beta));
            if (abs(cotarg) < 0.01 ) % Going to be singular
                %fprintf('Almost singular mode');
                eps = pi +  EDpm(plus, betaarg);
                if eps == 0
                    eps = 0.001;
                end
                if (abs(eps) >= 1)
                    msg = 'eps more than 1';
                    error(msg)
                end
                ret = n * exp(1i * pi / 4) * (sqrt(2 * pi * k * Lambda) * sign(eps) - 2 * k * Lambda * eps * exp(1i * pi / 4));
            else
                a = 2 * (cos(betaarg / 2) ^ 2);
                X = k * Lambda * a;
                if (X < 0.8)
                    F = sqrt(pi * X) * (1 - sqrt(X) / (0.7 * sqrt(X) + 1.2));
                else
                    F = 1 - 0.8 / ((X + 1.25) ^ 2);
                end
                F = F * (exp(1i * (pi / 4) * (1 - sqrt(X) / (X + 1.4))));
                ret = F / tan(cotarg);
            end
        end
        ret = equationquarter(beta, true) + equationquarter(beta, false);
    end
D = equationhalf(thetad - thetai) - EDpm(softboundary, equationhalf(thetad + thetai));
D = frontfactor * D;
end
%% Create biquad coefficients from ZPK parameters. Traceable by dlgradient

function [b, a] = BiquadCoefficients(zR, zI, pR, pI, k, numBiquads, numObservations)

    isDlarray = isdlarray(zR);

    if isDlarray
        [b, a] = deal(dlarray(zeros(3,numBiquads,numObservations)));
    else
        [b, a] = deal(zeros(3,numBiquads,numObservations));
    end
    for i = 1:numBiquads
        b(:,i,:) = [ones(1,numObservations); -2 * zR(i,:); zR(i,:).^2 + zI(i,:).^2];
        a(:,i,:) = [ones(1,numObservations); -2 * pR(i,:); pR(i,:).^2 + pI(i,:).^2];
    end
    b(:,1,:) = squeeze(b(:,1,:)) .* k;
    b = real(b);
    a = real(a);
    if isDlarray
        b = stripdims(b);
        a = stripdims(a);
    end

end
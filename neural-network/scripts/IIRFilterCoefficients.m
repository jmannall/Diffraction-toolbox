%% Create IIR filter coefficients from ZPK parameters. Traceable by dlgradient

function [b, a] = IIRFilterCoefficients(z, p, k, numIIRFilters, numObservations)

    isDlarray = isdlarray(z);

    if isDlarray
        [b, a] = deal(dlarray(zeros(2,numIIRFilters,numObservations)));
    else
        [b, a] = deal(zeros(2,numIIRFilters,numObservations));
    end
    for i = 1:numIIRFilters
        b(:,i,:) = [ones(1,numObservations); -z(i,:)];
        a(:,i,:) = [ones(1,numObservations); -p(i,:)];
    end
    b(:,1,:) = squeeze(b(:,1,:)) .* k;
    b = real(b);
    a = real(a);
    if isDlarray
        b = stripdims(b);
        a = stripdims(a);
    end
end
%% Add phase shift based on path length to transfer function
% Optional phase shift of pi determined by idx

function tfcomplex = PhaseShift(tfcomplex, pathLength, fvec, c, idx)
    
    delay = pathLength / c;
    mag = abs(tfcomplex);
    phase = angle(tfcomplex);

    phaseShift = -2 * pi * delay * fvec;
    phase = phase + phaseShift';
    tfcomplex = mag.*exp(phase*1i);

    if nargin > 4
        phase = angle(tfcomplex);
        phase(:,idx) = phase(:,idx) + pi;
        tfcomplex = mag.*exp(phase*1i);
    end
end
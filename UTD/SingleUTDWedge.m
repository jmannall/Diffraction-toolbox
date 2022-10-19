function [tfmag, fvec, tfcomplex] = SingleUTDWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phi, always, controlparameters, c)

    v = pi / deg2rad(wedgeIndex);
    thetaS = deg2rad(thetaS);
    thetaR = deg2rad(thetaR);
    
    nfft = controlparameters.nfft;
    fvec = controlparameters.fs / nfft * (0:nfft / 2 - 1);
    tfcomplex = zeros(1, nfft / 2);
    for x = 1:length(fvec)
        xm = EDutd_core(thetaS, radiusS, thetaR, radiusR, v, phi, 2 * pi * fvec(x) / c, false);
        if (~always)
            xm = xm * 2;
        end
        tfcomplex(x) = xm;
    end
    tfcomplex = tfcomplex';
    tfmag = mag2db(abs(tfcomplex));
end
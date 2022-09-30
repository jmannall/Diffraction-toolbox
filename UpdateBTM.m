function UpdateBTM(S)
    [y, ~] = S.fcn(S);   % @(S) BTMInf(S);
    figure(S.fh);
    SetString(S.wLabel,'Wedge',S.w);
    SetString(S.bALabel,'Bending Angle',S.bA);
    SetString(S.mALabel,'Min Angle',S.mA);
    set(S.LN, 'YData', y);
end
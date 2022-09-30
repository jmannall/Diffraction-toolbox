function S = CreateFirstOrderIIRSliders(S)

    % Inputs
    numZeros = length(S.z);
    numPoles = length(S.p);
    width = 500 / numPoles;
    textwidth = width;
    height = 40;
    spacing = 60;
    step = [0.01 0.01];
    callbackFcn = @SliderIIRCB;

    %% Z Sliders
    for i = 1:numZeros
        S.zLabel(i) = CreateLabel(10+spacing+(i-1)*(width+textwidth)+width, 120, textwidth, height);
        S.zSlider(i) = CreateSlider(spacing+(i-1)*(width+textwidth), 120, width, height, -1, 1, S.z(i), step, callbackFcn, 'z', i);
    end

    %% P Sliders
    for i = 1:numPoles
        S.pLabel(i) = CreateLabel(10+spacing+(i-1)*(width+textwidth)+width, 40, textwidth, height);
        S.pSlider(i) = CreateSlider(spacing+(i-1)*(width+textwidth), 40, width, height, -1, 1, S.p(i), step, callbackFcn, 'p', i);
    end

    %% K slider
    S.kLabel = CreateLabel(1710, 65, 300, 50);
    S.kSlider = CreateSlider(1400, 80, 300, 40, 0, 1, S.k, step, callbackFcn, 'k');

    %% Labels
    x = 1250;
    width = 500;
    height = 40;
    S.dcLabel = CreateLabel(x, 600, width, height);
    S.nyqLabel = CreateLabel(x, 520, width, height);
    S.shelfLabel = CreateLabel(x, 440, width, height);
    S.fcLLabel = CreateLabel(x, 360, width, height);
    S.fcHLabel = CreateLabel(x, 280, width, height);
end
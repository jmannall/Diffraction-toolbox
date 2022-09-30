function S = CreateBTMSliders(S)

    % Inputs
    width = 500;
    textwidth = width;
    height = 40;
    spacingX = 1200;
    wY = 800;
    bAY = 700;
    mAY = 600;
    step = [1/360 1/72];
    callbackFcn = @SliderBTMCB;

    %% Wedge index slider
    S.wLabel = CreateLabel(10+spacingX, 10+wY+height, textwidth, height);
    S.wSlider = CreateSlider(spacingX, wY, width, height, 0.001, 359.999, S.w, step, callbackFcn, 'w');

    %% Bending angle slider
    S.bALabel = CreateLabel(10+spacingX, 10+bAY+height, textwidth, height);
    S.bASlider = CreateSlider(spacingX, bAY, width, height, 0.001, 359.998, S.bA, step, callbackFcn, 'bA');

    %% Minimum angle slider
    S.mALabel = CreateLabel(10+spacingX, 10+mAY+height, textwidth, height);
    S.mASlider = CreateSlider(spacingX, mAY, width, height, 0.001, 179.999, S.mA, step, callbackFcn, 'mA');
end
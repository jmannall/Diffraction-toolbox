%% Create slider for an interctive figure

function slider = CreateSlider(x, y, width, height, min, max, value, step, callbackFcn, idx, i)

    if nargin < 11
        i = 1;
    end

    slider = uicontrol('style', 'slide', ...
        'unit', 'pixels', 'position', [x y width height], ...
        'min', min, 'max', max, 'value', value, ...
        'sliderstep', step, ...
        'callback', {callbackFcn, idx, i});
end
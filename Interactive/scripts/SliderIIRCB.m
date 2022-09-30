%% Callback function to update IIR sliders

function SliderIIRCB(mSlider, EventData, idx, i)
    S = guidata(mSlider);  % Get S struct from the figure
    S.(idx)(i) = get(mSlider, 'Value');
    S.update(S);
    guidata(mSlider, S);  % Store modified S in figure
end
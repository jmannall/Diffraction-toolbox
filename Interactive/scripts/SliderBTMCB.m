%% Callback function to update BTM sliders

function SliderBTMCB(mSlider, EventData, idx, i)
    S = guidata(mSlider);  % Get S struct from the figure
    newValue = max(get(mSlider, 'Min'), min(get(mSlider, 'Max'), round(get(mSlider, 'Value'))));
    set(mSlider, 'Value', newValue);
    S.(idx)(i) = newValue;
    S.update(S);
    guidata(mSlider, S);  % Store modified S in figure
end
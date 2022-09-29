function SliderCB_BTM(mSlider, EventData, Param1)
    % Callback for both sliders
    S = guidata(mSlider);  % Get S struct from the figure
    newValue = max(get(mSlider, 'Min'), min(get(mSlider, 'Max'), round(get(mSlider, 'Value'))));
    set(mSlider, 'Value', newValue);
    S.(Param1)= newValue;
    S.update(S);
    guidata(mSlider, S);  % Store modified S in figure
end
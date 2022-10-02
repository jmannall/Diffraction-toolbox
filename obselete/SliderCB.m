function SliderIIRCB(mSlider, EventData, Param1, Param2)
    % Callback for both sliders
    S = guidata(mSlider);  % Get S struct from the figure
    S.(Param1)(Param2) = get(mSlider, 'Value');
    S.update(S);
    guidata(mSlider, S);  % Store modified S in figure
end
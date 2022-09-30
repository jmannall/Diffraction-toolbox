function S = CreateIIRSliders(S,numZeros)
    width = 500 / numZeros;
    textwidth = width;
    height = 40;
    spacing = 60;
    for i = 1:numZeros
        S.zLabel(i) = uicontrol('Style','text','Position',[10+spacing+(i-1)*(width+textwidth)+width 120 textwidth height],'fontsize',20,'HorizontalAlignment','left');
        S.zSlider(i) = uicontrol('style','slide',...
                     'unit','pixels', 'position',[spacing+(i-1)*(width+textwidth) 120 width height],...
                     'min',-1,'max',1,'value', S.z(i),...
                     'sliderstep',[0.01 0.01],...
                     'callback', {@SliderCB, 'z', i});
    end
end
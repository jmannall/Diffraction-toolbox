function S = CreateP(S,numPoles)
    width = 500 / numPoles;
    textwidth = width;
    height = 40;
    spacing = 60;
    for i = 1:numPoles
        S.pLabel(i) = uicontrol('Style','text','Position',[10+spacing+(i-1)*(width+textwidth)+width 40 textwidth height],'fontsize',20,'HorizontalAlignment','left');
        S.pSlider(i) = uicontrol('style','slide',...
                     'unit','pixels', 'position',[spacing+(i-1)*(width+textwidth) 40 width height],...
                     'min',-1,'max',1,'value', S.p(i),...
                     'sliderstep',[0.01 0.01],...
                     'callback', {@SliderCB, 'p', i});
    end
end
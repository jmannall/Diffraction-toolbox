function S = CreateK(S)
    S.kLabel = uicontrol('Style','text','Position',[1710 65 300 50],'fontsize',20,'HorizontalAlignment','left');
    S.kSlider = uicontrol('style','slide',...
                 'unit','pixels', 'position',[1400 80 300 40],...
                 'min',0,'max',1,'value', S.k,...
                 'sliderstep',[0.01 0.01],...
                 'callback', {@SliderCB, 'k', 1});
end
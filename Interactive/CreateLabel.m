function label = CreateLabel(x, y, width, height)
    label = uicontrol('Style', 'text', 'Position', [x y width height], ...
        'FontSize', 20, 'HorizontalAlignment', 'left');
end
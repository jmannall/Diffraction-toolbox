function [lineOne, lineTwo] = CreateAnimatedLinePlot()

    figure
    C = colororder;
    lineOne = animatedline('Color',C(1,:));
    lineTwo = animatedline('Color',C(2,:));
    ylim([0 inf])
    xlabel("Iteration")
    ylabel("Loss")
    grid on
end